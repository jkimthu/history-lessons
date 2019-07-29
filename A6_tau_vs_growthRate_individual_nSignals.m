%% Figure A6: what nutrient signals are associated with each individual-level bin?


%  Goal: could nutrient signal explain the spread observed between
%        individuals in fluctuating environments?



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell interdivision time and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line




%  Last edit: jen, 2019 July 26
%  Commit: add mean, median and cv calculations

%  OK let's go!

%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. create array of experiments to use in analysis
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
compiled_tau = cell(length(exptArray),1);
compiled_birthLength = cell(length(exptArray),1);
compiled_lambda = cell(length(exptArray),1);
compiled_signal = cell(length(exptArray),1);


%% Part 2. collect single cell interdivision times and instantaneous growth rates

%  Strategy:
%
%       1.  loop through each experiment to collect data
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. for each condition in experiment...
%                       5. isolate condition specific data
%                       6. isolate curveDuration, timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. trim condition and growth rate data to include only full cell cycles
%                       9. isolate isDrop, timestamp, and timeSinceBirth data
%                      10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
%                      11. remove zeros, which occur if no full track data exists at a drop
%                      12. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      13. truncate data to stabilized regions
%                      14. if no div data in steady-state, skip condition
%                          else, trim outliers (those 3 std dev away from median) from final dataset
%                      15. bin growth rates and signals by cell cycle, to match organization of birth size data
%                      16. store condition data into one variable per experiment
%               17. store experiment data into single variable for further analysis
%      18. save hard earned data


% 1. loop through each experiment to collect data
for e = 1:length(exptArray)
 
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. initialize experiment specific variables
    interDiv = cell(length(bubbletime),1);
    mu_instantaneous = cell(length(bubbletime),1);
    nSignal = cell(length(bubbletime),1);
    
    
    
    % 4. load measured experiment data    
    %experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    %cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    

    
    % for each condition in experiment
    for condition = 1:length(bubbletime)
            
            
        % 5. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end date
        
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes isDrop trackNum timestamps_sec
        
        
        
        % 8. trim condition and growth rate data to include only full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        growthRates_fullOnly = growthRates(curveFinder > 0,:);
        curveIDs_fullOnly = curveFinder(curveFinder > 0,:);   % for trimming growth rate
        curveIDs_unique = unique(curveIDs_fullOnly);          % for assigning birth sizes
        clear conditionData curveFinder growthRates growthRates_all
        
        
        
        % 9. calculate binary nutrient signals
        [binaryNutrientSignal, nScore] = nutrientScore(timescale,conditionData_fullOnly);
        
        
        
        % 10. isolate timestamp, isDrop, and interdiv time data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        volume = getGrowthParameter(conditionData_fullOnly,'volume');     % calculated va_vals (cubic um)
        tau = getGrowthParameter(conditionData_fullOnly,'curveDurations'); % tau repeated for entire growth curve
        clear timestamps
        
        
        
        % 11. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthSize = volume(isDrop==1);
        final_tau = tau(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly timestamps_hr
        
        
        
        % 12. remove zeros, which occur if no full track data exists at a drop
        interdivisionT = final_tau(final_birthSize > 0);
        birthTimestamps = finalTimestamps(final_birthSize > 0);
        curveIDs = curveIDs_unique(final_birthSize > 0);
        clear final_birthSize finalTimestamps volume final_tau tau
        
        
        
        % 13. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            interdivisionT_bubbleTrimmed = interdivisionT(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            interdivisionT_bubbleTrimmed = interdivisionT;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear maxTime isDrop interdivisionT curveIDs birthTimestamps
        
        
        
        % 14. truncate data to stabilized regions
        minTime = 3;
        interdivisionT_fullyTrimmed = interdivisionT_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear interdivisionT_bubbleTrimmed  curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
        % 15. if no div data in steady-state, skip condition
        if isempty(interdivisionT_fullyTrimmed) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            tau_median = median(interdivisionT_fullyTrimmed);
            tau_std_temp = std(interdivisionT_fullyTrimmed);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            tau_temp = interdivisionT_fullyTrimmed(interdivisionT_fullyTrimmed <= (tau_median+tau_std_temp*3)); % cut largest vals, over 3 std out
            IDs_temp = curveIDs_fullyTrimmed_cc(interdivisionT_fullyTrimmed <= (tau_median+tau_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            tau_final = tau_temp(tau_temp >= (tau_median-tau_std_temp*3));          % cut smallest vals, over 3 std out
            IDs_final = IDs_temp(tau_temp >= (tau_median-tau_std_temp*3));   
            clear tau_median tau_std_temp tau_temp IDs_temp
            
            % iv. remove corresponding growth rates and signals from datasets
            trimmedIDs = setdiff(curveIDs_unique,IDs_final);    % curve IDs in growth rate dataset, NOT in final IDs trimmed by cell cycle
            toTrim = ismember(curveIDs_fullOnly,trimmedIDs);   % vector of what to trim or not in growth rate
            trimmed_curves_insta = curveIDs_fullOnly(toTrim == 0);
            trimmed_mus = growthRates_fullOnly(toTrim == 0);
            trimmed_signals = binaryNutrientSignal(toTrim == 0);
            clear toTrim trimmedIDs curveIDs_fullOnly growthRates_fullOnly binaryNutrientSignal
            
            
                 
            % 16. bin growth rates and nutrient signal by cell cycle, to match organization of birth size data
            mus_binned = accumarray(trimmed_curves_insta,trimmed_mus,[],@(x) {x});
            signals_binned = accumarray(trimmed_curves_insta,trimmed_signals,[],@(x) {x});
            mus = mus_binned(~cellfun('isempty',mus_binned));
            signals = signals_binned(~cellfun('isempty',signals_binned));
            clear trimmed_curves_insta trimmed_mus trimmed_signals
            
            
            
            % 16. remove cell cycles with negative mean growth rate
            meanGR = cellfun(@nanmean,mus);
            
            mus = mus(meanGR > 0);
            tau_final = tau_final(meanGR > 0);
            signals = signals(meanGR > 0);
            
            
            % 17. store condition data into one variable per experiment
            interDiv{condition} = tau_final;
            mu_instantaneous{condition} = mus;
            nSignal{condition} = signals;
        
            
        end
        
    end
      
    
    % 18. store experiment data into single variable for further analysis
    compiled_tau{e} = interDiv;
    compiled_lambda{e} = mu_instantaneous;
    compiled_signal{e} = nSignal;
end

% 19. save hard earned data
save('A6_data.mat','compiled_tau','compiled_lambda','compiled_signal','exptArray')


%% Part 3. bin cell cycles in 2D by tau and mean growth rate (lambda)

% goal: heated scatter where color represents nutrient signal!

% strategy: 
%
%   0. initialize complete meta data
%   0. initialize plotting parameters
%   0. initialize plotting parameters
%   1. loop through experiments to format data and plot
%   2. initialize experiment meta data
%   3. isolate experiment data
%   4. calculate lambda and nScores from compiled mus
%   5. remove data with too large growth rates
%   6. plot unheated scatter to ensure that heated scatter plots accurately
%   7. assign each cc to bins based on tau and growth rate
%   8. sort cell cycles into 2D bins: tau vs growth rate
%   9. plot scatter, heated by cell cycle counts
%  10. plot scatter, heated by mean nScore
%  11. plot scatter, heated by signal type
%  12. save plots



clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A6_data.mat')

% 0. initialize plotting parameters
%palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};
shape = 'o';

plot_tau = 90;
plot_lambda = 4;


% 0. initialize binning
binSize_mu = 0.25; % 1/h
binSize_tau = 1;  % min

max_tau = 110;
max_lambda = 5;


% 1. loop through experiments to format data and plot
for ee = 11:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(ee);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    for condition = 1%:4
        
        % 3. isolate experiment data
        eeTaus = compiled_tau{ee}{condition}./60; % convert from sec to h
        eeMus = compiled_lambda{ee}{condition};
        eeSignals = compiled_signal{ee}{condition};
        
        
        % 4. calculate lambda and nScores from compiled mus
        if isempty(eeTaus) == 1
            
            cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
            
            figure(1)
            plotName = strcat('A6-scatter-fig1-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(2)
            plotName = strcat('A6-scatter-fig2-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(3)
            plotName = strcat('A6-heated-cellCounts-fig3-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(4)
            plotName = strcat('A6-heated-cellCounts-fig4-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(5)
            plotName = strcat('A6-heated-nScore-fig5-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(6)
            plotName = strcat('A6-heated-nScore-fig6-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            continue
        else
            
            eeLambdas = cellfun(@nanmean,eeMus);
            eeNscores = cellfun(@nanmean,eeSignals);
            clear eeMus
            
            
            % 5. remove data with too large growth rates
            taus = eeTaus(eeLambdas <= max_lambda);
            nSignals = eeSignals(eeLambdas <= max_lambda);
            nScores = eeNscores(eeLambdas <= max_lambda);
            lambdas = eeLambdas(eeLambdas <= max_lambda);
            
            
            
            % 6. plot unheated scatter to ensure that heated scatter plots accurately
            color = rgb(palette(condition));
            
            % birth volume vs. mean growth rate
            figure(1)
            subplot(2,2,condition)
            plot(eeLambdas,eeTaus,'o','Color',color)
            hold on
            title(date)
            xlabel('mean growth rate (1/h)')
            ylabel('birth volume (cubic um)')
            axis([0 max_lambda 0 max_tau])
            
            clear color eeLambdas eeTaus eeSignals eeNscores
            
            
            
            % 7. assign each cc to bins based on volume/length and growth rate
            bin_tau = ceil(taus/binSize_tau);
            bin_lambda = ceil(lambdas/binSize_mu);
            
            
            % 8. sort cell cycles into 2D bins: volume/length vs growth rate
            
            %    i. initialize 2D bin matrix
            bin_matrix_tau = zeros(max_tau/binSize_tau,max_lambda/binSize_mu);           % vol vs lambda
            
            
            %   ii. assign cell cycle into a 2D bin, identified as a linear index of bin matrix
            linidx_tau = sub2ind(size(bin_matrix_tau),bin_tau,bin_lambda);

            
            %  iii. identify bins in matrix with cell cycle data, count number of cc per bin
            bins_unique_tau = unique(linidx_tau);
            out = [bins_unique_tau,histc(linidx_tau(:),bins_unique_tau)]; % outputs linear index in column 1, number of counts per index in column 2
            idx = out(:,1);
            counts = out(:,2);
            
            
            %  iv. assign count values to matrix bins with 2D subscripts
            bin_matrix_tau(out(:,1)) = out(:,2);
            [i,j] = ind2sub(size(bin_matrix_tau),idx);
            clear bin_matrix_tau out idx
            
            %   v. scale subscripts by bin size such to reflect real size and growth rate values
            binned_lambdas = j.*binSize_mu;
            binned_tau = i.*binSize_tau;
            clear j i
            
            
            % 9. plot scatter, heated by cell cycle counts
            figure(2)
            subplot(2,2,condition)
            scatter(binned_lambdas,binned_tau,60,counts,'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 plot_lambda 10 plot_tau])
            title(strcat(num2str(condition),'-',date))
            xlabel('mean growth rate (1/h)')
            ylabel('interdivison time (min)')
            
            
            
            
            % 10. plot scatter, heated by mean nScore
            nScores_by_bins = accumarray(linidx_tau,nScores,[],@(x) {x});
            nScores = nScores_by_bins(~cellfun('isempty',nScores_by_bins));
            nScores = cellfun(@mean,nScores);
            clear nScores_by_bins 
            
            
            figure(3)
            subplot(2,2,condition)
            scatter(binned_lambdas,binned_tau,60,nScores,'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 plot_lambda 10 plot_tau])
            title(strcat(num2str(condition),'-',date,'-heated-nScore'))
            xlabel('mean growth rate (1/h)')
            ylabel('interdiv time (min)')
            
            clear nScores
            
            
            % 11. plot scatter, heated by signal type
            %     i. determine types
            if timescale == 3600 && condition == 1
                
                signalType = zeros(length(nSignals),1);
                onlyH = 1;
                onlyL = 2;
                H2L = 3;
                L2H = 4;
                HLH = 5;
                LHL = 6;
                HLHL =7;
                LHLH = 8;
                
                for cc = 1:length(nSignals)
                    
                    currSignal = nSignals{cc};
                    ds_dt = diff(currSignal);
                    currShifts = length(find(ds_dt ~= 0));
                    
                    if currShifts == 0
                        if currSignal(1) == 1
                            signalType(cc) = onlyH;
                        else
                            signalType(cc) = onlyL;
                        end
                    elseif currShifts == 1
                        if currSignal(1) == 1
                            signalType(cc) = H2L;
                        else
                            signalType(cc) = L2H;
                        end
                    elseif currShifts == 2
                        if currSignal(1) == 1
                            signalType(cc) = HLH;
                        else
                            signalType(cc) = LHL;
                        end
                    elseif currShifts == 3
                        if currSignal(1) == 1
                            signalType(cc) = HLHL;
                        else
                            signalType(cc) = LHLH;
                        end
                    end
                end
                clear currSignal ds_dt currShifts cc
                
                types_by_bins = accumarray(linidx_tau,signalType,[],@(x) {x});
                types = types_by_bins(~cellfun('isempty',types_by_bins));
          
                
                
                % count signals of interest per bin
                sigArray = [onlyH; onlyL; H2L; L2H; HLH; LHL; HLHL; LHLH];
                sigNames = {'onlyH','onlyL','H2L','L2H','HLH','LHL','HLHL','LHLH'};
                signal_counts = zeros(length(types),length(sigArray));

                
                for sg = 1:length(sigArray)
                    
                    currSig = sigArray(sg);
                    currName = sigNames{sg};
                    
                    % binned by volume data
                    for ii = 1:length(types)
                        currBin = types{ii};
                        numSig = find(currBin == currSig);
                        if isempty(numSig) == 1
                            signal_counts(ii,sg) = 0;
                        else
                            signal_counts(ii,sg) = length(numSig);
                        end
                    end
                    clear numSig currBin ii
                    
                end
                clear currSig sg
                
                
                % calculate fraction of signals of interest
                counts_8div = [counts,counts,counts,counts,counts,counts,counts,counts];
                signal_frac_v = signal_counts./counts_8div;
                clear counts_8div
                
                
                for sgg = 1:length(sigArray)
                    
                    % plot fraction of current signal of interest
                    figure(4)
                    subplot(2,2,condition)
                    scatter(binned_lambdas,binned_tau,60,signal_frac_v(:,sgg),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
                    colormap(parula) % see colormap documentation
                    colorbar;
                    axis([0 plot_lambda 10 plot_tau])
                    title(strcat(num2str(condition),'-',date,'-',sigNames{sgg}))
                    xlabel('mean growth rate (1/h)')
                    ylabel('tau (cubic um)')
                    
                    
                    % save plots for only condition 1
                    figure(4)
                    plotName = strcat('A6-heated-',sigNames{sgg},'-',date,'-vol');
                    saveas(gcf,plotName,'epsc')
                    close(gcf)

                    
                    % for sanity check, see below
                    types_counted(sgg) = length(find(signalType == sigArray(sgg)));
                    
                end
                
                % sanity check, plot bar plot of signal types
                figure(5)
                bar(sigArray(1:sgg),types_counted)
                ylabel('counts')
                xlabel('types')
                title(date)
                set(gca,'xticklabel',sigNames)
                
            end
            
            
        end
        
        
        % 12. save plots
%         cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
%         
%         figure(1)
%         plotName = strcat('A6-scatter-fig1-',date,'-tau');
%         saveas(gcf,plotName,'epsc')
%         close(gcf)
%         
%         figure(2)
%         plotName = strcat('A6-heated-cellCounts-fig2-',date,'-tau');
%         saveas(gcf,plotName,'epsc')
%         close(gcf)
%         
%         figure(3)
%         plotName = strcat('A6-heated-nScore-fig3-',date,'-tau');
%         saveas(gcf,plotName,'epsc')
%         close(gcf)
%         
        
    end
    types_counted
    
    if timescale == 3600
        
%         figure(5)
%         plotName = strcat('A6-bar-fig5-',date);
%         saveas(gcf,plotName,'epsc')
%         close(gcf)
%         
    end
    
end


%% Part 4. plot distribution of interdivision time, comparing 15 and 60 min conditions


clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A6_data.mat')

% 0. initialize plotting parameters
%palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};


% 1. isolate conditions of interest, convert all from sec to min

devs_rep1 = compiled_tau{11}{1}./60; % 60 min, 2019-01-29
devs_rep2 = compiled_tau{12}{1}./60; % 60 min, 2019-01-31
devs_rep3 = compiled_tau{13}{1}./60; % 60 min, 2019-02-01

nots_rep1 = compiled_tau{8}{1}./60; % 15 min, 2017-01-12
nots_rep2 = compiled_tau{9}{1}./60; % 15 min, 2018-01-16
nots_rep3 = compiled_tau{10}{1}./60; % 15 min, 2018-01-17
nots_rep4 = compiled_tau{7}{1}./60; % 15 min, 2017-11-13


% 2. remove un-physiological data: interdivision times < 10 min
devs_rep1 = devs_rep1(devs_rep1 > 9);
devs_rep2 = devs_rep2(devs_rep2 > 9);
devs_rep3 = devs_rep3(devs_rep3 > 9);

nots_rep1 = nots_rep1(nots_rep1 > 9);
nots_rep2 = nots_rep2(nots_rep2 > 9);
nots_rep3 = nots_rep3(nots_rep3 > 9);
nots_rep4 = nots_rep4(nots_rep4 > 9);


% 3. compile into cell array of "left" and "right" bins
deviants = {devs_rep1, devs_rep2, devs_rep3};
goodies = {nots_rep1, nots_rep2, nots_rep3};
goodies_all = {nots_rep4, nots_rep1, nots_rep2, nots_rep3};

% 4. plot replicate distributions of interdivision time as violin plot
figure(1)
distributionPlot(deviants,'widthDiv',[2 1],'histOri','left','color',rgb('DarkBlue'),'showMM',2)
hold on
distributionPlot(gca,goodies,'widthDiv',[2 2],'histOri','right','color',rgb('DeepSkyBlue'),'showMM',2)
xlabel('replicate')
ylabel('inter-division time (min)')
title('histograms of interdivision time across conditions')
legend('left: 60 min','right: 15 min')


% 5. report mean and median of each replicate
deviants_mean = cellfun(@mean,deviants);
deviants_median = cellfun(@median,deviants);
deviants_std = cellfun(@std,deviants);
deviants_cv = deviants_std./deviants_mean * 100

goodies_mean = cellfun(@mean,goodies_all);
goodies_median = cellfun(@median,goodies_all);
goodies_std = cellfun(@std,goodies_all);
goodies_cv = goodies_std./goodies_mean * 100







