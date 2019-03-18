%% C2: traits binned and averaged by nutrient phase at birth


% Goal: are growth phenotypes predicted by when a cell is born?
%       plot various traits after binning birth data by nutrient phase
%
%           1. interdivision time
%           2. mean growth rate across cell cycle
%           3. volume at birth
%           4. volume at division
%           5. volume added
%           6. division size divided by birth size
%           7. nScore
%           8. cell cycle at time of shift

%       repeat with subpopulations: 1 std above and and below division size


%  Strategy:
%
%     0. initialize complete meta data
%     0. for experiments of interest...
%








%  Last edit: Jen Nguyen, 2019 Mar 12

%  Commit: first commit,


% ok let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = [13;14;15]; % use corresponding dataIndex values


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
palette_below = {'LightSkyBlue','BlueViolet','Gold','LightCoral'};


% 0. initialize binning parameters
%binsPerPeriod = 12;
timePerBin = 2.5; % min 



%%

% 1. for all experiments in dataset
for e = 1:length(exptArray)
    
    
    
    % 2. collect experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    

    % 2b. calculate bins per period
    binsPerPeriod = (timescale/60)/timePerBin;
    
    
    % 3. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    % 4. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType filename
    
    
    
    % for each condition in experiment...
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        color_above = rgb(palette(condition));
        color_below = rgb(palette_below(condition));
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % col 11 = calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % col 2  = timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % col 4  = isDrop, 1 marks a birth event
        curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % col 20 = track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveID,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes timestamps_sec isDrop trackNum
        
        
        
        % 8. trim condition and growth rate data to include only full cell cycles
        fullData = conditionData(curveID > 0,:);
        growthRates_fullOnly = growthRates(curveID > 0,:);
        clear conditionData curveID growthRates
        
        
        
        % 9. isolate volume, isDrop, curveID and timestamp data
        volumes = getGrowthParameter(fullData,'volume');
        isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
        curveIDs = getGrowthParameter(fullData,'curveFinder');
        timestamps_sec = getGrowthParameter(fullData,'timestamp');
        timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
        clear timestamps_sec
        
        
        % 10. identity unique cell cycles by ID number
        cellCycles = curveIDs(isDrop == 1);
        birthTimes = timestamps_hr(isDrop == 1);
        clear fullData
        
        
        % 11. remove birth times prior to 3 hr
        birthTimes_post3 = birthTimes(birthTimes > 3);
        cellCycles_post3 = cellCycles(birthTimes > 3);
        clear cellCycles birthTimes
        
        
        % 12. remove birth times post bubbles
        birthTimes_final = birthTimes_post3(birthTimes_post3 < bubbletime(condition));
        cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
        clear cellCycles_post3 birthTimes_post3 birthTimes_final
        
        
        % 13. for remaining cell cycles, identify:
        %       1. volume at birth
        %       2. volume at division
        %       3. interdivision time
        %       4. mean growth rate
        %       5. nutrient phase at birth
        ccData = nan(length(cellCyles_final),9);
        
        for cc = 1:length(cellCyles_final)
            
            % isolate data for current cell cycle
            currentVolumes = volumes(curveIDs == cellCyles_final(cc));
            currentTimestamps = timestamps_hr(curveIDs == cellCyles_final(cc));
            currentGrowthRates = growthRates_fullOnly(curveIDs == cellCyles_final(cc));
            
            ccData(cc,1) = currentVolumes(1);       % V_birth
            ccData(cc,2) = currentVolumes(end);     % V_division
            ccData(cc,3) = currentVolumes(end) - currentVolumes(1); % added volume
            ccData(cc,4) = currentTimestamps(1);    % T_birth
            ccData(cc,5) = currentTimestamps(end);  % T_division
            ccData(cc,6) = currentTimestamps(end) - currentTimestamps(1); % interdivision time
            ccData(cc,7) = nanmean(currentGrowthRates); % mean growth rate
            ccData(cc,8) = nanstd(currentGrowthRates);  % standard deviation in growth rate
            
            
            % for nutrient phase
            % i.  re-define period to begin at start of low nutrient pulse, by
            %     subtracting quarter period from corrected timestamp
            currentTimestamps_shifted = (currentTimestamps * 3600) - (timescale/4); % converted from hr to sec
            
            % ii. assign elements of timestamp vector to period fraction bins
            timeInPeriods = currentTimestamps_shifted/timescale; % unit = sec/sec
            timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
            assignedBin = ceil(timeInPeriodFraction * binsPerPeriod);
            
            ccData(cc,9) = assignedBin(1); % nutrient phase at birth (in bin)
            
        end
        clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly
        clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final
        clear assignedBin timeInPeriodFraction timeInPeriods currentTimestamps_shifted
        
        
        
        % 14. trim cell cycle data to avoid 1 & 2 point cell cycles
        addedVol = ccData(:,3);
        intermediate = ccData(addedVol > 0,:);
        
        grDevs = intermediate(:,8);
        traits = intermediate(grDevs > 0,:);
        clear ccData addedVol grDevs intermediate
        
        
        
        % 15. exclude outliers based on cell size (both birthsize and divsize)
        %     a) id median and std of division size
        %     b) id median and std of birth size
        %     c) find indeces in both vectors that are within 3 std
        %     d) keep data from indeces in both vectors
        V_birth = traits(:,1);
        V_div = traits(:,2);
        
        divSize_median = median(V_div);
        divSize_std = std(V_div);
        
        birthSize_median = median(V_birth);
        birthSize_std = std(V_birth);
        
        div_bigOutlier = find(V_div > (divSize_median+divSize_std*3));
        div_smallOutlier = find(V_div < (divSize_median-divSize_std*3));
        div_outliers = [div_bigOutlier; div_smallOutlier];
        
        birth_bigOutlier = find(V_birth > (birthSize_median+birthSize_std*3));
        birth_smallOutlier = find(V_birth < (birthSize_median-birthSize_std*3));
        birth_outliers = [birth_bigOutlier; birth_smallOutlier];
        
        V_division_binary = ones(length(V_div),1);
        V_division_binary(div_outliers) = 0;
        
        V_birth_binary = ones(length(V_birth),1);
        V_birth_binary(birth_outliers) = 0;
        
        V_summed = V_division_binary + V_birth_binary;
        
        traits_final = traits(V_summed == 2,:);
        
        clear birthSize_median birthSize_std divSize_median divSize_std
        clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
        clear V_birth_binary V_division_binary birth_outliers div_outliers V_div V_birth V_summed interdiv
        
        
        
        
        % 16. bin data and normalize bin counts by total counts

        % some conditions do not have data after trimming
        if isempty(traits_final) == 1
            
            continue
            
        else
            
            % A. split condition data into two populations: above and below average
            divSize = traits_final(:,2);
            divSize_mean = mean(divSize);
            divSize_std = std(divSize);
            
            sizeThresh_above = divSize_mean + divSize_std;
            sizeThres_below = divSize_mean - divSize_std;
            
            traits_above = traits_final(divSize > sizeThresh_above,:);
            traits_below = traits_final(divSize < sizeThres_below,:);
            clear traits divSize divSize_mean divSize_std sizeThresh_above sizeThres_below
            
            
            % B. bins traits of interest from subpop "above" by birth phase
            binNum = traits_above(:,9);
            binNum_b = traits_below(:,9);

            % above
            tau = accumarray(binNum, traits_above(:,6), [], @(x) {x}); %   1. interdivision time
            mu = accumarray(binNum, traits_above(:,7), [], @(x) {x});  %   2. mean growth rate across cell cycle
            birthSize = accumarray(binNum, traits_above(:,1), [], @(x) {x}); %   3. volume at birth
            divSize = accumarray(binNum, traits_above(:,2), [], @(x) {x}); %   4. volume at division
            addedSize = accumarray(binNum, traits_above(:,3), [], @(x) {x});  %   5. volume added
            sizeRatio = accumarray(binNum, traits_above(:,2)./traits_above(:,1), [], @(x) {x}); %   6. division size divided by birth size
            
            % below
            tau_b = accumarray(binNum_b, traits_below(:,6), [], @(x) {x}); %   1. interdivision time
            mu_b = accumarray(binNum_b, traits_below(:,7), [], @(x) {x});  %   2. mean growth rate across cell cycle
            birthSize_b = accumarray(binNum_b, traits_below(:,1), [], @(x) {x}); %   3. volume at birth
            divSize_b = accumarray(binNum_b, traits_below(:,2), [], @(x) {x}); %   4. volume at division
            addedSize_b = accumarray(binNum_b, traits_below(:,3), [], @(x) {x});  %   5. volume added
            sizeRatio_b = accumarray(binNum_b, traits_below(:,2)./traits_below(:,1), [], @(x) {x}); %   6. division size divided by birth size
            
            % if shorter than full period, extend data vector
            if length(tau_b) < binsPerPeriod
                
                tau_b{binsPerPeriod} = [];
                mu_b{binsPerPeriod} = [];
                birthSize_b{binsPerPeriod} = [];
                divSize_b{binsPerPeriod} = [];
                addedSize_b{binsPerPeriod} = [];
                sizeRatio_b{binsPerPeriod} = [];
                
            end
            
            binCounts{1,e}{condition} = cellfun(@length,tau);
            binCounts{2,e}{condition} = cellfun(@length,tau_b);
            
          
            
            % C. plot mean and s.e.m. of cell cycle duration over nutrient period
            %    i. convert bin # to absolute time (sec)
            periodInBins = linspace(1, binsPerPeriod, binsPerPeriod);
            periodInTime = timePerBin*periodInBins';
            clear periodInBins
            
            
            % 11. (ii) repeat quarter period on both sides and plot over period fraction
%             quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
%             quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
%             quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
%             
%             signal_quarterOne = durations_binned_mean(quarterOne,1);
%             signal_quarterFour = durations_binned_mean(quarterFour,1);
%             
%             time_quarterOne = periodInTime(quarterOne,1) + timescale;
%             time_quarterZero = quarterZero*timePerBin;
            
%             mean_durationSignal_stitched = [signal_quarterFour; durations_binned_mean; signal_quarterOne];
%             time_stitched = [time_quarterZero'; periodInTime; time_quarterOne];
            
            figure(1) % interdivision time (tau, above)
            hold on
            errorbar(periodInTime, cellfun(@mean,tau)*60 , cellfun(@std,tau)*60, 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,tau_b)*60 , cellfun(@std,tau_b)*60, 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('interdivision time (min)')
            title('tau vs. nutrient phase of birth')
            axis([0,65,0,80])

            
            figure(2) % mean growth rate (mu, above)
            hold on
            errorbar(periodInTime, cellfun(@mean,mu), cellfun(@std,mu), 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,mu_b), cellfun(@std,mu_b), 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('growth rate (1/h)')
            title('mu vs. nutrient phase of birth')
            axis([0,65,0,4.5])
            
            
            figure(3) % birth size
            errorbar(periodInTime, cellfun(@mean,birthSize), cellfun(@std,birthSize), 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,birthSize_b), cellfun(@std,birthSize_b), 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('size at birth (um^3)')
            title('birth size vs. nutrient phase of birth')
            axis([0,65,0,8])
            

            figure(4) % division size
            errorbar(periodInTime, cellfun(@mean,divSize), cellfun(@std,divSize), 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,divSize_b), cellfun(@std,divSize_b), 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('size at division (um^3)')
            title('division size vs. nutrient phase of birth')
            axis([0,65,2,15])
            
            
            figure(5) % added size
            errorbar(periodInTime, cellfun(@mean,addedSize), cellfun(@std,addedSize), 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,addedSize_b), cellfun(@std,addedSize_b), 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('added size (um^3)')
            title('added size vs. nutrient phase of birth')
            axis([0,65,0,10])
            
            
            figure(6) % size ratio
            errorbar(periodInTime, cellfun(@mean,sizeRatio), cellfun(@std,sizeRatio), 'Marker','o','Color',color_above)
            hold on
            errorbar(periodInTime, cellfun(@mean,sizeRatio_b), cellfun(@std,sizeRatio_b), 'Marker','o','Color',color_below)
            xlabel('time (min)')
            ylabel('Vdiv/Vbirth')
            title('Volume ratio vs. nutrient phase of birth')
            axis([0,65,0,4])
            
            
            
        end
        
        
        
    end
    
%     % 19. save plots
    cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
    
    figure(1)
    plotName = strcat('C2-fig1-tau-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    plotName = strcat('C2-fig2-mu-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(3)
    plotName = strcat('C2-fig3-vBirth-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(4)
    plotName = strcat('C2-fig4-vDiv-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(5)
    plotName = strcat('C2-fig5-addedV-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(6)
    plotName = strcat('C2-fig6-vRatio-v-nPhase-',date,'-timePerBin-2p5');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    clc
    
    
end

%% figure 7: birth events per nutrient phase bin

% sum births below 1 st dev
a = binCounts{2,1}{1};
b = binCounts{2,2}{1};
c = binCounts{2,3}{1};
d = a + b + c;

% sun births above 1 st dev
e = binCounts{1,1}{1};
f = binCounts{1,2}{1};
g = binCounts{1,3}{1};
h = e + f + g;

figure(7)
stem(d)
hold on
stem(h,'Color',[ 0.9100 0.4100 0.1700])
xlabel('bin no.')
ylabel('birth events')
legend('below','above')
title('births over nutrient period')


% sum births in low
clear a b c d e f g h
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
palette_below = {'LightSkyBlue','BlueViolet','Gold','LightCoral'};

for s = 2:4
    
    color = rgb(palette(s));
    color_b = rgb(palette_below(s));
    
    s1 = binCounts{1,1}{s};
    s2 = binCounts{1,2}{s};
    s3 = binCounts{1,3}{s};
    steady_sum = s1 + s2 + s3;
    
    s4 = binCounts{2,1}{s};
    s5 = binCounts{2,2}{s};
    s6 = binCounts{2,3}{s};
    steady_sum_below = s4 + s5 + s6;
    
    figure(s)
    stem(steady_sum_below,'Color',color_b)
    hold on
    stem(steady_sum,'Color',color)
    xlabel('bin no.')
    ylabel('birth events')
    legend('below','above')
    title('births over nutrient period')
    ylim([0 210])

    
    
end
