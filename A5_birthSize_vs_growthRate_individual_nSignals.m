%% Figure A5: what nutrient signals are associated with each individual-level bin?


%  Goal: could nutrient signal explain the spread observed between
%        individuals in fluctuating environments?



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell birth size and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line




%  Last edit: jen, 2019 July 22
%  Commit: first commit, throw nutrient signal into the mix


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
compiled_birthVolume = cell(length(exptArray),1);
compiled_birthLength = cell(length(exptArray),1);
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
%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
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
    birthSize = cell(length(bubbletime),1);
    birthLength = cell(length(bubbletime),1);
    mu_instantaneous = cell(length(bubbletime),1);
    nSignal = cell(length(bubbletime),1);
    
    
    
    % 4. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
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
        %clear timescale
        
        
        
        % 10. isolate timestamp, isDrop, length and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        volumes = getGrowthParameter(conditionData_fullOnly,'volume');     % calculated va_vals (cubic um)
        majorAxis = getGrowthParameter(conditionData_fullOnly,'length');   % length (um)
        clear timestamps
        
        
        
        % 11. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthSize = volumes(isDrop==1);
        final_birthLength = majorAxis(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly
        
        
        
        % 12. remove zeros, which occur if no full track data exists at a drop
        Vbirth = final_birthSize(final_birthSize > 0);
        Lbirth = final_birthLength(final_birthSize > 0);
        birthTimestamps = finalTimestamps(final_birthSize > 0);
        curveIDs = curveIDs_unique(final_birthSize > 0);
        clear final_birthSize finalTimestamps volumes majorAxis
        
        
        
        % 13. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            Vbirth_bubbleTrimmed = Vbirth(birthTimestamps <= maxTime,:);
            Lbirth_bubbleTrimmed = Lbirth(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            Vbirth_bubbleTrimmed = Vbirth;
            Lbirth_bubbleTrimmed = Lbirth;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear timestamps_hr maxTime isDrop Vbirth Lbirth curveIDs birthTimestamps
        
        
        
        % 14. truncate data to stabilized regions
        minTime = 3;
        Vbirth_fullyTrimmed = Vbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        Lbirth_fullyTrimmed = Lbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear Vbirth_bubbleTrimmed Lbirth_bubbleTrimmed curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
        % 15. if no div data in steady-state, skip condition
        if isempty(Vbirth_fullyTrimmed) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            birthSize_median = median(Vbirth_fullyTrimmed);
            birthSize_std_temp = std(Vbirth_fullyTrimmed);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            birthSize_temp = Vbirth_fullyTrimmed(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3)); % cut largest vals, over 3 std out
            birthLength_temp = Lbirth_fullyTrimmed(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3));
            IDs_temp = curveIDs_fullyTrimmed_cc(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            birthSize_final = birthSize_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));          % cut smallest vals, over 3 std out
            birthLength_final = birthLength_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3)); 
            IDs_final = IDs_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));   
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
            
            
            
            % 16. remove cell cycles with negative mean growth rate
            meanGR = cellfun(@nanmean,mus);
            
            mus = mus(meanGR > 0);
            birthSize_final = birthSize_final(meanGR > 0);
            birthLength_final = birthLength_final(meanGR > 0);
            signals = signals(meanGR > 0);
            
            
            % 17. store condition data into one variable per experiment
            birthSize{condition} = birthSize_final;
            birthLength{condition} = birthLength_final;
            mu_instantaneous{condition} = mus;
            nSignal{condition} = signals;
        
            
        end
        
    end
      
    
    % 18. store experiment data into single variable for further analysis
    compiled_birthVolume{e} = birthSize;
    compiled_birthLength{e} = birthLength;
    compiled_lambda{e} = mu_instantaneous;
    compiled_signal{e} = nSignal;
end

% 19. save hard earned data
save('A5_data.mat','compiled_birthVolume','compiled_birthLength','compiled_lambda','compiled_signal','exptArray')


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
load('A5_data.mat')

% 0. initialize plotting parameters
%palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
palette = {'Indigo','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};
shape = 'o';


% 0. initialize binning
binSize_mu = 0.25; % 1/h
binSize_vol = 0.25;  % cubic um
binSize_length = 0.25; % um

max_vol = 15;
max_length = 11;
max_lambda = 5;


% 1. loop through experiments to format data and plot
for ee = 9:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(ee);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    for condition = 1:4
        
        % 3. isolate experiment data
        eeVols = compiled_birthVolume{ee}{condition};
        eeLengths = compiled_birthLength{ee}{condition};
        eeMus = compiled_lambda{ee}{condition};
        eeSignals = compiled_signal{ee}{condition};
        
        
        % 4. calculate lambda and nScores from compiled mus
        eeLambdas = cellfun(@nanmean,eeMus);
        eeNscores = cellfun(@nanmean,eeSignals);
        clear eeMus
        
        
        % 5. remove data with too large growth rates
        volumes = eeVols(eeLambdas <= max_lambda);
        lengths = eeLengths(eeLambdas <= max_lambda);
        nSignals = eeSignals(eeLambdas <= max_lambda);
        nScores = eeNscores(eeLambdas <= max_lambda);
        lambdas = eeLambdas(eeLambdas <= max_lambda);
        
        
        
        % 6. plot unheated scatter to ensure that heated scatter plots accurately
        color = rgb(palette(condition));
        
        % birth volume vs. mean growth rate
        figure(1)
        subplot(2,2,condition)
        plot(eeLambdas,eeVols,'o','Color',color)
        hold on
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        axis([0 max_lambda 0 max_vol])
        
        % birth length vs. mean growth rate
        figure(2)
        subplot(2,2,condition)
        plot(eeLambdas,eeLengths,'o','Color',color)
        hold on
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        axis([0 max_lambda 0 max_length])
        clear color
        clear eeLambdas eeVols eeLengths eeSignals eeNscores
        
        
        
        % 7. assign each cc to bins based on volume/length and growth rate
        bin_birthVol = ceil(volumes/binSize_vol);
        bin_birthLength = ceil(lengths/binSize_length);
        bin_lambda = ceil(lambdas/binSize_mu);
        
        
        % 8. sort cell cycles into 2D bins: volume/length vs growth rate
        
        %    i. initialize 2D bin matrix
        bin_matrix_vol = zeros(max_vol/binSize_vol,max_lambda/binSize_mu);           % vol vs lambda
        bin_matrix_length = zeros(max_length/binSize_length,max_lambda/binSize_mu);  % length vs lambda
        
        %   ii. assign cell cycle into a 2D bin, identified as a linear index of bin matrix
        linidx_vol = sub2ind(size(bin_matrix_vol),bin_birthVol,bin_lambda);
        linidx_len = sub2ind(size(bin_matrix_length),bin_birthLength,bin_lambda);
        
        %  iii. identify bins in matrix with cell cycle data, count number of cc per bin
        bins_unique_vol = unique(linidx_vol);
        out_v = [bins_unique_vol,histc(linidx_vol(:),bins_unique_vol)]; % outputs linear index in column 1, number of counts per index in column 2
        idx_v = out_v(:,1); 
        counts_v = out_v(:,2);
        
        bins_unique_len = unique(linidx_len);
        out_l = [bins_unique_len,histc(linidx_len(:),bins_unique_len)];
        idx_l = out_l(:,1);
        counts_l = out_l(:,2);
        
        %  iv. assign count values to matrix bins with 2D subscripts
        bin_matrix_vol(out_v(:,1)) = out_v(:,2);
        bin_matrix_length(out_l(:,1)) = out_l(:,2);
        
        [iv,jv] = ind2sub(size(bin_matrix_vol),idx_v);
        [il,jl] = ind2sub(size(bin_matrix_length),idx_l);
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
        axis([0 max_lambda 0 max_vol])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        
        figure(4)
        subplot(2,2,condition)
        scatter(binned_lambdas_l,binned_length,60,counts_l,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 max_lambda 0 max_length])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        

        
        % 10. plot scatter, heated by mean nScore
        nScores_by_vbins = accumarray(linidx_vol,nScores,[],@(x) {x});
        nScores_by_lbins = accumarray(linidx_len,nScores,[],@(x) {x});
        
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
        axis([0 max_lambda 0 max_vol])
        title(strcat(num2str(condition),'-',date,'-heated-nScore'))
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        
        figure(6)
        subplot(2,2,condition)
        scatter(binned_lambdas_l,binned_length,60,nScores_len,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 max_lambda 0 max_length])
        title(strcat(num2str(condition),'-',date))
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        
        clear nScores_len nScores_vol
        
        
        % 11. plot scatter, heated by signal type
        %     i. determine types
        if timescale == 3600
            
            signalType = zeros(length(nSignals),1);
            onlyH = 1;
            onlyL = 2;
            H2L = 3;
            L2H = 4;
            HLH = 5;
            LHL = 6;
            
            for cc = 1:length(nSignals)
                
                currSignal = nSignals{cc};
                ds_dt = diff(currSignal);
                currShifts = length(find(ds_dt ~= 0));
                currShift_types = ds_dt(ds_dt ~= 0);
                
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
                end
            end
            clear currSignal ds_dt currShifts currShift_types cc
            
            types_by_vbins = accumarray(linidx_vol,signalType,[],@(x) {x});
            types_by_lbins = accumarray(linidx_len,signalType,[],@(x) {x});
            
            types_vol = types_by_vbins(~cellfun('isempty',types_by_vbins));
            types_len = types_by_lbins(~cellfun('isempty',types_by_lbins));
            
           
            % count signals of interest per bin
            sigArray = [H2L; L2H; onlyH; onlyL];
            signal_counts_v = zeros(length(types_vol),length(sigArray));
            signal_counts_l = zeros(length(types_len),length(sigArray));
            
            for sg = 1:length(sigArray)
                
                currSig = sigArray(sg);
                
                % binned by volume data
                for ii = 1:length(types_vol)
                    currBin = types_vol{ii};
                    numSig = find(currBin == currSig);
                    if isempty(numSig) == 1
                        signal_counts_v(ii,currSig) = 0;
                    else
                        signal_counts_v(ii,currSig) = length(numSig);
                    end
                end
                
                % binned by length data
                for ii = 1:length(types_len)
                    currBin = types_len{ii};
                    numSig = find(currBin == currSig);
                    if isempty(numSig) == 1
                        signal_counts_l(ii,currSig) = 0;
                    else
                        signal_counts_l(ii,currSig) = length(numSig);
                    end
                end
                 
            end
            clear currSig numSig currBin ii
            
            % calculate fraction of signals of interest
            counts_v_4div = [counts_v,counts_v,counts_v,counts_v];
            counts_l_4div = [counts_l,counts_l,counts_l,counts_l];
            signal_frac_v = signal_counts_v./counts_v_4div;
            signal_frac_l = signal_counts_l./counts_l_4div;
            clear counts_v_4div counts_l_4div
            
            
            % plot fraction of H2L
            figure(7)
            subplot(2,2,condition)
            scatter(binned_lambdas_v,binned_vol,60,signal_frac_v(:,1),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_vol])
            title(strcat(num2str(condition),'-',date,'-fracH2L'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth volume (cubic um)')
            
            figure(8)
            subplot(2,2,condition)
            scatter(binned_lambdas_l,binned_length,60,signal_frac_l(:,1),'filled')
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_length])
            title(strcat(num2str(condition),'-',date,'-fracH2L'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth length (um)')
            
            
            % plot fraction of L2H
            figure(9)
            subplot(2,2,condition)
            scatter(binned_lambdas_v,binned_vol,60,signal_frac_v(:,2),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_vol])
            title(strcat(num2str(condition),'-',date,'-fracL2H'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth volume (cubic um)')
            
            figure(10)
            subplot(2,2,condition)
            scatter(binned_lambdas_l,binned_length,60,signal_frac_l(:,2),'filled')
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_length])
            title(strcat(num2str(condition),'-',date,'-fracL2H'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth length (um)')
            
            % plot fraction of only H
            figure(11)
            subplot(2,2,condition)
            scatter(binned_lambdas_v,binned_vol,60,signal_frac_v(:,3),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_vol])
            title(strcat(num2str(condition),'-',date,'-fracOnlyH'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth volume (cubic um)')
            
            figure(12)
            subplot(2,2,condition)
            scatter(binned_lambdas_l,binned_length,60,signal_frac_l(:,3),'filled')
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_length])
            title(strcat(num2str(condition),'-',date,'-fraconlyH'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth length (um)')
            
            
            % plot fraction of only L
            figure(13)
            subplot(2,2,condition)
            scatter(binned_lambdas_v,binned_vol,60,signal_frac_v(:,4),'filled') % scatter(x axis bin, y axis bin, circle size, value of color)
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_vol])
            title(strcat(num2str(condition),'-',date,'-fracOnlyL'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth volume (cubic um)')
            
            figure(14)
            subplot(2,2,condition)
            scatter(binned_lambdas_l,binned_length,60,signal_frac_l(:,4),'filled')
            colormap(parula) % see colormap documentation
            colorbar;
            axis([0 max_lambda 0 max_length])
            title(strcat(num2str(condition),'-',date,'-fraconlyL'))
            xlabel('mean growth rate (1/h)')
            ylabel('birth length (um)')
            
            % sanity check, plot bar plot of signal types
            types_all = [onlyH, onlyL, H2L, L2H, HLH, LHL];
            types_counted(onlyH) = length(find(signalType == onlyH));
            types_counted(onlyL) = length(find(signalType == onlyL));
            types_counted(H2L) = length(find(signalType == H2L));
            types_counted(L2H) = length(find(signalType == L2H));
            types_counted(HLH) = length(find(signalType == HLH));
            types_counted(LHL) = length(find(signalType == LHL));
            
            figure(15)
            bar(types_all,types_counted)
            ylabel('counts')
            xlabel('types')
            title(date)
            barNames = {'onlyH'; 'onlyL'; 'H2L'; 'L2H'; 'HLH';'LHL' };
            set(gca,'xticklabel',barNames)
        end
        
        
        
    end
    
    
     % 12. save plots
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
    
    if timescale == 3600
        figure(7)
        plotName = strcat('A5-heated-H2L-fig7-',date,'-vol');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(8)
        plotName = strcat('A5-heated-H2L-fig8-',date,'-length');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(9)
        plotName = strcat('A5-heated-L2H-fig9-',date,'-vol');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(10)
        plotName = strcat('A5-heated-L2H-fig10-',date,'-length');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(11)
        plotName = strcat('A5-heated-onlyH-fig11-',date,'-vol');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(12)
        plotName = strcat('A5-heated-onlyH-fig12-',date,'-length');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(13)
        plotName = strcat('A5-heated-onlyL-fig13-',date,'-vol');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(14)
        plotName = strcat('A5-heated-onlyL-fig14-',date,'-length');
        saveas(gcf,plotName,'epsc')
        close(gcf)
        
        figure(15)
        plotName = strcat('A5-bar-fig15-',date);
        saveas(gcf,plotName,'epsc')
        close(gcf)
    end
    
end







