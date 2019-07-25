%% Figure A7: what nutrient signals are associated with each individual-level bin?


%  Goal: could nutrient signal explain the spread observed between
%        individuals in fluctuating environments?



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell birth size and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line




%  Last edit: jen, 2019 July 25
%  Commit: first commit to look at added size as a function of nutrient signal


%  OK let's go!

%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. create array of experiments to use in analysis
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
compiled_addedVolume = cell(length(exptArray),1);
compiled_addedLength = cell(length(exptArray),1);
compiled_lambda = cell(length(exptArray),1);
compiled_signal = cell(length(exptArray),1);


%% Part 2. collect single cell birth size and instantaneous growth rates

%  Strategy:
%
%       1.  loop through each experiment to collect data
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. for each condition in experiment...
%                       5. isolate condition specific data
%                       6. isolate added volume, added length, timestamp, drop, curve, and trackNum data
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
for e = 11:length(exptArray)
 
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. initialize experiment specific variables
    addedVolume = cell(length(bubbletime),1);
    addedLength = cell(length(bubbletime),1);
    mu_instantaneous = cell(length(bubbletime),1);
    nSignal = cell(length(bubbletime),1);
    
    
    
    % 4. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    

    
    % for each condition in experiment
    for condition = 1%:length(bubbletime)
            
            
        % 5. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end date
        
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        del_V = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates_all = calculateGrowthRate(del_V,timestamps_sec,isDrop,curveFinder,trackNum);
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
        
        
        
        % 10. isolate timestamp, isDrop, length and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');           % isDrop, 1 marks a birth event
        addedV = getGrowthParameter(conditionData_fullOnly,'addedVA');         % added volume repeated across full cycle
        addedL = getGrowthParameter(conditionData_fullOnly,'addedLength');   % added length (um)
        clear timestamps
        
        
        
        % 11. extract only one value for added size, use birth event
        final_addedV = addedV(isDrop==1);
        final_addedL = addedL(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly
        
        
        
        % 12. remove zeros, which occur if no full track data exists at a drop
        delta_v = final_addedV(final_addedV > 0);
        delta_l = final_addedL(final_addedV > 0);
        birthTimestamps = finalTimestamps(final_addedV > 0);
        curveIDs = curveIDs_unique(final_addedV > 0);
        clear final_addedV final_addedL finalTimestamps addedV addedL
        
        
        
        % 13. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            deltaV_bubbleTrimmed = delta_v(birthTimestamps <= maxTime,:);
            deltaL_bubbleTrimmed = delta_l(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            deltaV_bubbleTrimmed = delta_v;
            deltaL_bubbleTrimmed = delta_l;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear timestamps_hr maxTime isDrop delta_v delta_l curveIDs birthTimestamps
        
        
        
        % 14. truncate data to stabilized regions
        minTime = 3;
        deltaV_fullyTrimmed = deltaV_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        deltaL_fullyTrimmed = deltaL_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear deltaV_bubbleTrimmed deltaL_bubbleTrimmed curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
        % 15. if no div data in steady-state, skip condition
        if isempty(deltaV_fullyTrimmed) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            dV_median = median(deltaV_fullyTrimmed);
            dV_std_temp = std(deltaV_fullyTrimmed);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            dV_temp = deltaV_fullyTrimmed(deltaV_fullyTrimmed <= (dV_median+dV_std_temp*3)); % cut largest vals, over 3 std out
            dL_temp = deltaL_fullyTrimmed(deltaV_fullyTrimmed <= (dV_median+dV_std_temp*3));
            IDs_temp = curveIDs_fullyTrimmed_cc(deltaV_fullyTrimmed <= (dV_median+dV_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            dV_final = dV_temp(dV_temp >= (dV_median-dV_std_temp*3));          % cut smallest vals, over 3 std out
            dL_final = dL_temp(dV_temp >= (dV_median-dV_std_temp*3)); 
            IDs_final = IDs_temp(dV_temp >= (dV_median-dV_std_temp*3));   
            clear birthSize_median birthSize_std_temp birthSize_temp birthLength_temp IDs_temp
            
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
            
            
            
            % 17. remove cell cycles with negative mean growth rate
            meanGR = cellfun(@nanmean,mus);
            
            mus = mus(meanGR > 0);
            dV_final = dV_final(meanGR > 0);
            dL_final = dL_final(meanGR > 0);
            signals = signals(meanGR > 0);
            
            
            % 18. store condition data into one variable per experiment
            addedVolume{condition} = dV_final;
            addedLength{condition} = dL_final;
            mu_instantaneous{condition} = mus;
            nSignal{condition} = signals;
        
            
        end
        
    end
      
    
    % 18. store experiment data into single variable for further analysis
    compiled_addedVolume{e} = addedVolume;
    compiled_addedLength{e} = addedLength;
    compiled_lambda{e} = mu_instantaneous;
    compiled_signal{e} = nSignal;
end

% 19. save hard earned data
save('A7_data.mat','compiled_addedVolume','compiled_addedLength','compiled_lambda','compiled_signal','exptArray')


%% Part 3. bin cell cycles in 2D by birth size and mean growth rate (lambda)

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
%   7. assign each cc to bins based on volume/length and growth rate
%   8. sort cell cycles into 2D bins: volume/length vs growth rate
%   9. plot scatter, heated by cell cycle counts
%  10. plot scatter, heated by mean nScore
%  11. plot scatter, heated by signal type
%  12. save plots



clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A7_data.mat')

% 0. initialize plotting parameters
%palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};
shape = 'o';

plot_length = 6;
plot_vol = 8;
plot_lambda = 4;


% 0. initialize binning
binSize_mu = 0.25; % 1/h
binSize_vol = 0.25;  % cubic um
binSize_length = 0.25; % um

max_vol = 15;
max_length = 11;
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
        eeAV = compiled_addedVolume{ee}{condition};
        eeAL = compiled_addedLength{ee}{condition};
        eeMus = compiled_lambda{ee}{condition};
        eeSignals = compiled_signal{ee}{condition};
        
        
        % 4. calculate lambda and nScores from compiled mus
        if isempty(eeAV) == 1
            
            cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
            
            figure(1)
            plotName = strcat('A5-scatter-fig1-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(2)
            plotName = strcat('A5-scatter-fig2-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(3)
            plotName = strcat('A5-heated-cellCounts-fig3-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(4)
            plotName = strcat('A5-heated-cellCounts-fig4-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(5)
            plotName = strcat('A5-heated-nScore-fig5-',date,'-vol');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            figure(6)
            plotName = strcat('A5-heated-nScore-fig6-',date,'-length');
            saveas(gcf,plotName,'epsc')
            close(gcf)
            
            continue
        end
        eeLambdas = cellfun(@nanmean,eeMus);
        eeNscores = cellfun(@nanmean,eeSignals);
        clear eeMus
        
        
        % 5. remove data with too large growth rates
        del_V = eeAV(eeLambdas <= max_lambda);
        del_L = eeAL(eeLambdas <= max_lambda);
        nSignals = eeSignals(eeLambdas <= max_lambda);
        nScores = eeNscores(eeLambdas <= max_lambda);
        lambdas = eeLambdas(eeLambdas <= max_lambda);
        
        
        
        
        % 6. plot unheated scatter to ensure that heated scatter plots accurately
        color = rgb(palette(condition));
        
        % birth volume vs. mean growth rate
        figure(1)
        subplot(2,2,condition)
        plot(eeLambdas,eeAV,'o','Color',color)
        hold on
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        axis([0 max_lambda 0 max_vol])
        
        % birth length vs. mean growth rate
        figure(2)
        subplot(2,2,condition)
        plot(eeLambdas,eeAL,'o','Color',color)
        hold on
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        axis([0 max_lambda 0 max_length])
        clear color
        clear eeLambdas eeAV eeAL eeSignals eeNscores
        
        
        
        % 7. assign each cc to bins based on volume/length and growth rate
        bin_AV = ceil(del_V/binSize_vol);
        bin_AL = ceil(del_L/binSize_length);
        bin_lambda = ceil(lambdas/binSize_mu);
        
        
        
        % 8. remove added lengths than bin to zero
        bin_AV = bin_AV(bin_AL > 0);
        bin_AL = bin_AL(bin_AL > 0);
        bin_lambda = bin_lambda(bin_AL > 0);
        nScores = nScores(bin_AL > 0);
        nSignals = nSignals(bin_AL > 0);
        
        
        
        % 9. sort cell cycles into 2D bins: volume/length vs growth rate
        
        %    i. initialize 2D bin matrix
        bin_matrix_av = zeros(max_vol/binSize_vol,max_lambda/binSize_mu);           % vol vs lambda
        bin_matrix_al = zeros(max_length/binSize_length,max_lambda/binSize_mu);  % length vs lambda
        
        %   ii. assign cell cycle into a 2D bin, identified as a linear index of bin matrix
        linidx_av = sub2ind(size(bin_matrix_av),bin_AV,bin_lambda);
        linidx_al = sub2ind(size(bin_matrix_al),bin_AL,bin_lambda);
        
        %  iii. identify bins in matrix with cell cycle data, count number of cc per bin
        bins_unique_vol = unique(linidx_av);
        out_v = [bins_unique_vol,histc(linidx_av(:),bins_unique_vol)]; % outputs linear index in column 1, number of counts per index in column 2
        idx_v = out_v(:,1); 
        counts_v = out_v(:,2);
        
        bins_unique_len = unique(linidx_al);
        out_l = [bins_unique_len,histc(linidx_al(:),bins_unique_len)];
        idx_l = out_l(:,1);
        counts_l = out_l(:,2);
        
        %  iv. assign count values to matrix bins with 2D subscripts
        bin_matrix_av(out_v(:,1)) = out_v(:,2);
        bin_matrix_al(out_l(:,1)) = out_l(:,2);
        
        [iv,jv] = ind2sub(size(bin_matrix_av),idx_v);
        [il,jl] = ind2sub(size(bin_matrix_al),idx_l);
        clear bin_matrix_vol bin_matrix_length out_l idx_l out_v idx_v
        
        %   v. scale subscripts by bin size such to reflect real size and growth rate values
        binned_lambdas_v = jv.*binSize_mu;
        binned_vol = iv.*binSize_vol;
        
        binned_lambdas_l = jl.*binSize_mu;
        binned_length = il.*binSize_length;
        clear jv jl iv il
        
      
        % 9. plot scatter, heated by cell cycle counts
        figure(3)
        subplot(2,2,condition)
        scatter(binned_lambdas_v,binned_vol,60,counts_v,'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 plot_lambda 0 plot_vol])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('added volume (cubic um)')
        
        figure(4)
        subplot(2,2,condition)
        scatter(binned_lambdas_l,binned_length,60,counts_l,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 plot_lambda 0 plot_length])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('added length (um)')
        

        
        % 10. plot scatter, heated by mean nScore
        nScores_by_vbins = accumarray(linidx_av,nScores,[],@(x) {x});
        nScores_by_lbins = accumarray(linidx_al,nScores,[],@(x) {x});
        
        nScores_vol = nScores_by_vbins(~cellfun('isempty',nScores_by_vbins));
        nScores_len = nScores_by_lbins(~cellfun('isempty',nScores_by_lbins));
        
        nScores_vol = cellfun(@mean,nScores_vol);
        nScores_len = cellfun(@mean,nScores_len);
        clear nScores_by_vbins nScores_by_lbins
        
        
        figure(5)
        subplot(2,2,condition)
        scatter(binned_lambdas_v,binned_vol,60,nScores_vol,'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 plot_lambda 0 plot_vol])
        title(strcat(num2str(condition),'-',date,'-heated-nScore'))
        xlabel('mean growth rate (1/h)')
        ylabel('added volume (cubic um)')
        
        figure(6)
        subplot(2,2,condition)
        scatter(binned_lambdas_l,binned_length,60,nScores_len,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 plot_lambda 0 plot_length])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('added length (um)')
        
        clear nScores_len nScores_vol
        
        
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
            
            types_by_vbins = accumarray(linidx_av,signalType,[],@(x) {x});
            types_by_lbins = accumarray(linidx_al,signalType,[],@(x) {x});
            
            types_vol = types_by_vbins(~cellfun('isempty',types_by_vbins));
            types_len = types_by_lbins(~cellfun('isempty',types_by_lbins));
            
           
            % count signals of interest per bin
            sigArray = [onlyH; onlyL; H2L; L2H; HLH; LHL; HLHL; LHLH];
            sigNames = {'onlyH','onlyL','H2L','L2H','HLH','LHL','HLHL','LHLH'};
            signal_counts_v = zeros(length(types_vol),length(sigArray));
            signal_counts_l = zeros(length(types_len),length(sigArray));
            
            for sg = 1:length(sigArray)
                
                currSig = sigArray(sg);
                currName = sigNames{sg};
                
                % binned by volume data
                for ii = 1:length(types_vol)
                    currBin = types_vol{ii};
                    numSig = find(currBin == currSig);
                    if isempty(numSig) == 1
                        signal_counts_v(ii,sg) = 0;
                    else
                        signal_counts_v(ii,sg) = length(numSig);
                    end
                end
                
                % binned by length data
                for ii = 1:length(types_len)
                    currBin = types_len{ii};
                    numSig = find(currBin == currSig);
                    if isempty(numSig) == 1
                        signal_counts_l(ii,sg) = 0;
                    else
                        signal_counts_l(ii,sg) = length(numSig);
                    end
                end
                clear numSig currBin ii
                 
            end
            clear currSig sg
            
            % calculate fraction of signals of interest
            counts_v_8div = [counts_v,counts_v,counts_v,counts_v,counts_v,counts_v,counts_v,counts_v];
            counts_l_8div = [counts_l,counts_l,counts_l,counts_l,counts_l,counts_l,counts_l,counts_l];
            signal_frac_v = signal_counts_v./counts_v_8div;
            signal_frac_l = signal_counts_l./counts_l_8div;
            clear counts_v_8div counts_l_8div
            
            
            for sgg = 1:length(sigArray)
                
                % plot fraction of current signal of interest
                figure(7)
                subplot(2,2,condition)
                scatter(binned_lambdas_v,binned_vol,60,signal_frac_v(:,sgg),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
                colormap(parula) % see colormap documentation
                colorbar;
                axis([0 plot_lambda 0 plot_vol])
                title(strcat(num2str(condition),'-',date,'-',sigNames{sgg}))
                xlabel('mean growth rate (1/h)')
                ylabel('added volume (cubic um)')
                
                
                figure(8)
                subplot(2,2,condition)
                scatter(binned_lambdas_l,binned_length,60,signal_frac_l(:,sgg),'filled')
                colormap(parula) % see colormap documentation
                colorbar;
                axis([0 plot_lambda 0 plot_length])
                title(strcat(num2str(condition),'-',date,'-',sigNames{sgg}))
                xlabel('mean growth rate (1/h)')
                ylabel('added length (um)')
                
                
                % save plots for only condition 1
                figure(7)
                plotName = strcat('A7-heated-',sigNames{sgg},'-',date,'-av');
                saveas(gcf,plotName,'epsc')
                close(gcf)
                
                figure(8)
                plotName = strcat('A7-heated-',sigNames{sgg},'-',date,'-al');
                saveas(gcf,plotName,'epsc')
                close(gcf)
               
                
                % for sanity check, see below
                types_counted(sgg) = length(find(signalType == sigArray(sgg)));
                
            end
           
            % sanity check, plot bar plot of signal types
            figure(15)
            bar(sigArray(1:sgg),types_counted)
            ylabel('counts')
            xlabel('types')
            title(date)
            set(gca,'xticklabel',sigNames)
            
        end
        

    end
    
    
    % 12. save plots
    cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    
    figure(1)
    plotName = strcat('A7-scatter-fig1-',date,'-av');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    plotName = strcat('A7-scatter-fig2-',date,'-al');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(3)
    plotName = strcat('A7-heated-cellCounts-fig3-',date,'-av');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(4)
    plotName = strcat('A7-heated-cellCounts-fig4-',date,'-al');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(5)
    plotName = strcat('A7-heated-nScore-fig5-',date,'-av');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(6)
    plotName = strcat('A5-heated-nScore-fig6-',date,'-al');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    if timescale == 3600
 
        figure(15)
        plotName = strcat('A7-bar-fig15-',date);
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
    end
    
end







