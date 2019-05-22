%% Figure I: trait comparisons between steady-state and 0/100% nutrient signals 

%  Goal: for a given experiment,
%
%        0.  in fluctuating condition:
%               1.  compute nScore for all cell cycles after first 3h
%               2.  bin cells with 0% or 100% nutrient
%               3.  plot bar plot with mean and st dev of each trait
%        4.  in steady nutrient conditions:
%               5.  plot bar plot with mean and st dev of each trait, using
%                   data from all cell cycles after first 3h



%  Traits of interest: 
%
%       a) interdivision time
%       b) division size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score
%       f) added size




%  Last edit: jen, 2019 Apr 15
%  Commit: add mann-whitney u test between comparable bars


%  OK let's go!

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
exptArray = 13:15; % use corresponding dataIndex values



% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};



% 0. initialize nutrient signal classifications
classRules = [ 0,0,0,0,0; ... % 1  only low
    1,0,0,0,0; ...            % 2
    1,1,0,0,0; ...            % 3
    1,1,1,0,0; ...            % 4
    1,1,1,1,0; ...            % 5
    1,1,1,1,1; ...            % 6   only high
    0,1,1,1,1; ...            % 7
    0,0,1,1,1; ...            % 8
    0,0,0,1,1; ...            % 9
    0,0,0,0,1; ...            % 10
    1,0,0,0,1; ...            % 11
    1,0,0,1,1; ...            % 12
    1,1,0,0,1; ...            % 13
    1,1,0,1,1; ...            % 14  impossible?
    0,1,1,1,0; ...            % 15
    0,1,1,0,0; ...            % 16
    0,0,1,1,0; ...            % 17
    0,0,1,0,0; ...            % 18  impossible?
    1,1,1,0,1; ...            % 19
    1,0,1,1,1; ...            % 20
    0,1,0,0,0; ...            % 21
    0,0,0,1,0; ...            % 22
    0,1,1,0,1; ...            % 23
    1,0,1,0,1; ...            % 24
    0,1,0,1,0; ...            % 25
    0,1,0,0,1; ...            % 26
    0,1,0,1,1; ...            % 27
    1,0,0,1,0; ...            % 28
    0,0,1,0,1; ...            % 29
    1,0,0,1,1; ...            % 30
    ];


% 0. initialize array for concatenation of final cell cycles from each experimental dataset
traits_all = [];


%% Part 1. collect and concatenate cell cycle data from each fluctuating experiment

condition = 1;

% 1. for all experiments in dataset
for e = 1:length(exptArray)
    
    
    % 2. collect experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    timescale = storedMetaData{index}.timescale;
    xys = storedMetaData{index}.xys;
    disp(strcat(date, ': analyze!'))
    

    
    % 3. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    % 4. compile experiment data matrix
    xy_start = xys(condition,1);
    xy_end = xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType filename
    
    
    
    
    % 5. isolate volume (Va), timestamp, drop, curve, and trackNum data
    volumes = getGrowthParameter(conditionData,'volume');            % col 11 = calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % col 2  = timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');             % col 4  = isDrop, 1 marks a birth event
    curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');         % col 20 = track number (not ID from particle tracking)
    
    
    
    % 6. calculate growth rate
    growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveID,trackNum);
    growthRates = growthRates_all(:,specificColumn);
    clear volumes timestamps_sec isDrop trackNum
    
    
    
    % 7. trim condition and growth rate data to include only full cell cycles
    fullData = conditionData(curveID > 0,:);
    growthRates_fullOnly = growthRates(curveID > 0,:);
    
    
    
    % 8. calculate binary nutrient signals
    [binaryNutrientSignal, nScore] = nutrientScore(timescale,fullData);
    clear curveID growthRates growthRates_all
    
    
    
    % 9. isolate volume, isDrop, curveID, trackNum and timestamp data
    volumes = getGrowthParameter(fullData,'volume');
    isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
    curveIDs = getGrowthParameter(fullData,'curveFinder');
    trackNum = getGrowthParameter(fullData,'trackNum');
    timestamps_sec = getGrowthParameter(fullData,'timestamp');
    timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
    clear timestamps_sec
    
    
    % 10. identity unique cell cycles by ID number
    cellCycles = curveIDs(isDrop == 1);
    birthTimes = timestamps_hr(isDrop == 1);
    tracks = trackNum(isDrop == 1);
    clear fullData
    
    
    % 11. remove birth times prior to 3 hr
    birthTimes_post3 = birthTimes(birthTimes > 3);
    cellCycles_post3 = cellCycles(birthTimes > 3);
    tracks_post3 = tracks(birthTimes > 3);
    clear cellCycles birthTimes tracks
    
    
    % 12. remove birth times post bubbles
    cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
    tracks_final = tracks_post3(birthTimes_post3 < bubbletime(condition));
    clear cellCycles_post3 birthTimes_post3 tracks_post3
    
    
    % 13. for remaining cell cycles, identify:
    %       1. volume at birth
    %       2. volume at division
    %       3. interdivision time
    %       4. mean growth rate
    ccData = nan(length(cellCyles_final),11);
    
    for cc = 1:length(cellCyles_final)
        
        % isolate data for current cell cycle
        currentVolumes = volumes(curveIDs == cellCyles_final(cc));
        currentTimestamps = timestamps_hr(curveIDs == cellCyles_final(cc));
        currentGrowthRates = growthRates_fullOnly(curveIDs == cellCyles_final(cc));
        currentBinarySignal = binaryNutrientSignal(curveIDs == cellCyles_final(cc));
        currentNscore = nScore(curveIDs == cellCyles_final(cc));
        currentTrack = tracks_final(cc);
        currentCC = cellCyles_final(cc);
        
        if length(unique(currentNscore)) ~= 1
            error('Nscore in current cell cycle is not unique: error in calling cell cycle or score')
        end
        
        ccData(cc,1) = currentVolumes(1);       % V_birth
        ccData(cc,2) = currentVolumes(end);     % V_division
        ccData(cc,3) = currentVolumes(end) - currentVolumes(1); % added volume
        ccData(cc,4) = currentTimestamps(1);    % T_birth
        ccData(cc,5) = currentTimestamps(end);  % T_division
        ccData(cc,6) = currentTimestamps(end) - currentTimestamps(1); % interdivision time
        ccData(cc,7) = nanmean(currentGrowthRates); % mean growth rate
        ccData(cc,8) = nanstd(currentGrowthRates);  % standard deviation in growth rate
        ccData(cc,9) = currentNscore(1);               % nScore
        
        
        ccSignal{cc,1} = currentBinarySignal;        % binary signal (1 = high, 0 = low)
        
        
        isSwitch = [0; diff(currentBinarySignal)] ~= 0;
        numShifts = sum(isSwitch);
        if numShifts > 0
            switches{cc,1} = find(isSwitch == 1);
            shiftStage = switches{cc}(1)/length(currentBinarySignal); % fraction into cell cycle of first switch
        else
            switches{cc,1} = NaN;
            shiftStage = NaN;
        end
        
        
        ccData(cc,10) = numShifts;
        ccData(cc,11) = shiftStage; % first if multiple
        ccData(cc,12) = currentTrack;
        ccData(cc,13) = currentCC;
        
        
        % classify binary signal
        
        % bin cell cycle into 5ths,
        % considering 5ths as "high N" if over half includes high nutrient
        fractions = linspace(1,length(currentBinarySignal),length(currentBinarySignal))/length(currentBinarySignal);
        bins = ceil(fractions'*5);
        binned = round(accumarray(bins,currentBinarySignal,[],@mean));
        
        % classify signal
        for cl = 1:length(classRules)
            currentClass = classRules(cl,:);
            if isequal(binned',currentClass) == 1
                ccData(cc,14) = cl;
                break
            end
        end
        
        
    end
    clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly
    clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final
    clear binaryNutrientSignal nScore currentBinarySignal currentNscore
    clear isSwitch numShifts shiftStage currentTrack trackNum tracks_final currentCC
    clear bins fractions binned currentClass cl
    
    
    % 14. trim cell cycle data to avoid single point cell cycles
    addedVol = ccData(:,3);
    traits = ccData(addedVol > 0,:);
    signal = ccSignal(addedVol > 0,:);
    clear addedVol
    
    
    
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
    signal_final = signal(V_summed == 2,:);
    
    clear birthSize_median birthSize_std divSize_median divSize_std
    clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
    clear V_birth_binary V_division_binary birth_outliers div_outliers V_div V_birth V_summed interdiv
    clear switches
    
    
    
    % 16. exclude interdivision times that are less than 10 min long
    % note: still need to justify exclusion of interdiv times < 10 visually
    interdivs = traits_final(:,6)*60;
    traits_10plus = traits_final(interdivs > 10,:);
    clear signal_final traits_final signal traits interdivs
    
    
    % 17. concatenate data across replicates
    traits_all = [traits_all; traits_10plus];
    size(traits_all)
    
end

clear bubbletime ans traits_10plus xys date
clear ccData ccSignal conditionData experimentFolder


%% Part 2. identify 0%, 50% and 100% cell cycles by nScore

wiggle = 2.5; % percent error

% 1. bin traits by nScore
nScores = traits_all(:,9)*100;

traits_zero = traits_all(nScores < wiggle,:);
traits_100 = traits_all(nScores > (100 - wiggle),:);

traits_temp = traits_all(nScores > (50 - wiggle),:);
nScores_temp = nScores(nScores > (50 - wiggle),:);
traits_50 = traits_temp(nScores_temp < (50 + wiggle),:);
       
clear traits_all traits_temp nScores nScores_temp wiggle
clear classRules condition index e


%% Part 3. collect and concatenate cell cycle data from each steady condition/experiment

for condition = 2:4
    
    traits_cond = [];
    
    % 1. for all experiments in dataset
    for e = 1:length(exptArray)
        
        
        % 2. collect experiment meta data
        index = exptArray(e);
        date = storedMetaData{index}.date;
        expType = storedMetaData{index}.experimentType;
        bubbletime = storedMetaData{index}.bubbletime;
        timescale = storedMetaData{index}.timescale;
        xys = storedMetaData{index}.xys;
        disp(strcat(date, ': analyze!'))
        
        
        
        % 3. load measured experiment data
        experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
        cd(experimentFolder)
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
        load(filename,'D5','T');
        
        
        % 4. compile experiment data matrix
        xy_start = xys(condition,1);
        xy_end = xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear D5 T xy_start xy_end expType filename
        
        
        
        % 5. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % col 11 = calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % col 2  = timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % col 4  = isDrop, 1 marks a birth event
        curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % col 20 = track number (not ID from particle tracking)
        
        
        
        % 6. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveID,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes timestamps_sec isDrop trackNum
        
        
        
        % 7. trim condition and growth rate data to include only full cell cycles
        fullData = conditionData(curveID > 0,:);
        growthRates_fullOnly = growthRates(curveID > 0,:);
        
        
        
        % 8. assign binary nutrient signals
        if condition == 2
            binaryNutrientSignal = zeros(length(growthRates_fullOnly),1);
            nScore = zeros(length(growthRates_fullOnly),1);
        elseif condition == 3
            binaryNutrientSignal = ones(length(growthRates_fullOnly),1)*0.5;
            nScore = ones(length(growthRates_fullOnly),1)*0.5;
        else
            binaryNutrientSignal = ones(length(growthRates_fullOnly),1);
            nScore = ones(length(growthRates_fullOnly),1);
        end
        %[binaryNutrientSignal, nScore] = nutrientScore(timescale,fullData);
        clear curveID growthRates growthRates_all
        
        
        
        % 9. isolate volume, isDrop, curveID, trackNum and timestamp data
        volumes = getGrowthParameter(fullData,'volume');
        isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
        curveIDs = getGrowthParameter(fullData,'curveFinder');
        trackNum = getGrowthParameter(fullData,'trackNum');
        timestamps_sec = getGrowthParameter(fullData,'timestamp');
        timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
        clear timestamps_sec
        
        
        % 10. identity unique cell cycles by ID number
        cellCycles = curveIDs(isDrop == 1);
        birthTimes = timestamps_hr(isDrop == 1);
        tracks = trackNum(isDrop == 1);
        clear fullData
        
        
        % 11. remove birth times prior to 3 hr
        birthTimes_post3 = birthTimes(birthTimes > 3);
        cellCycles_post3 = cellCycles(birthTimes > 3);
        tracks_post3 = tracks(birthTimes > 3);
        clear cellCycles birthTimes tracks
        
        
        % 12. remove birth times post bubbles
        cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
        tracks_final = tracks_post3(birthTimes_post3 < bubbletime(condition));
        clear cellCycles_post3 birthTimes_post3 tracks_post3
        
        
        % 13. for remaining cell cycles, identify:
        %       1. volume at birth
        %       2. volume at division
        %       3. interdivision time
        %       4. mean growth rate
        ccData = nan(length(cellCyles_final),11);
        
        for cc = 1:length(cellCyles_final)
            
            % isolate data for current cell cycle
            currentVolumes = volumes(curveIDs == cellCyles_final(cc));
            currentTimestamps = timestamps_hr(curveIDs == cellCyles_final(cc));
            currentGrowthRates = growthRates_fullOnly(curveIDs == cellCyles_final(cc));
            currentBinarySignal = binaryNutrientSignal(curveIDs == cellCyles_final(cc));
            currentNscore = nScore(curveIDs == cellCyles_final(cc));
            currentTrack = tracks_final(cc);
            currentCC = cellCyles_final(cc);
            
            if length(unique(currentNscore)) ~= 1
                error('Nscore in current cell cycle is not unique: error in calling cell cycle or score')
            end
            
            ccData(cc,1) = currentVolumes(1);       % V_birth
            ccData(cc,2) = currentVolumes(end);     % V_division
            ccData(cc,3) = currentVolumes(end) - currentVolumes(1); % added volume
            ccData(cc,4) = currentTimestamps(1);    % T_birth
            ccData(cc,5) = currentTimestamps(end);  % T_division
            ccData(cc,6) = currentTimestamps(end) - currentTimestamps(1); % interdivision time
            ccData(cc,7) = nanmean(currentGrowthRates); % mean growth rate
            ccData(cc,8) = nanstd(currentGrowthRates);  % standard deviation in growth rate
            ccData(cc,9) = currentNscore(1);               % nScore
            
            
            ccSignal{cc,1} = currentBinarySignal;        % binary signal (1 = high, 0 = low)
            
            
            isSwitch = [0; diff(currentBinarySignal)] ~= 0;
            numShifts = sum(isSwitch);
            if numShifts > 0
                switches{cc,1} = find(isSwitch == 1);
                shiftStage = switches{cc}(1)/length(currentBinarySignal); % fraction into cell cycle of first switch
            else
                switches{cc,1} = NaN;
                shiftStage = NaN;
            end
            
            
            ccData(cc,10) = numShifts;
            ccData(cc,11) = shiftStage; % first if multiple
            ccData(cc,12) = currentTrack;
            ccData(cc,13) = currentCC;
            
            
            % classify binary signal as NaN (it's steady!)
            ccData(cc,14) = NaN;

            
        end
        clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly
        clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final
        clear binaryNutrientSignal nScore currentBinarySignal currentNscore
        clear isSwitch numShifts shiftStage currentTrack trackNum tracks_final currentCC
        clear bins fractions binned currentClass cl
        
        
        
        % 14. trim cell cycle data to avoid single point cell cycles
        addedVol = ccData(:,3);
        traits = ccData(addedVol > 0,:);
        signal = ccSignal(addedVol > 0,:);
        clear addedVol
        
        
        
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
        signal_final = signal(V_summed == 2,:);
        
        clear birthSize_median birthSize_std divSize_median divSize_std
        clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
        clear V_birth_binary V_division_binary birth_outliers div_outliers V_div V_birth V_summed interdiv
        clear switches
        
        
        
        % 16. exclude interdivision times that are less than 10 min long
        % note: still need to justify exclusion of interdiv times < 10 visually
        interdivs = traits_final(:,6)*60;
        traits_10plus = traits_final(interdivs > 10,:);
        clear signal_final traits_final signal traits interdivs
        
        
        % 17. concatenate data across replicates
        traits_cond = [traits_cond; traits_10plus];
        size(traits_cond)
        
    end
    

    % 18. store mean, st dev and counts of all cell cycles in condition
    traits_steady{condition}.values = traits_cond;
    traits_steady{condition}.means = mean(traits_cond);
    traits_steady{condition}.stds = std(traits_cond);
    traits_steady{condition}.count = length(traits_cond);
    

    
    clear bubbletime ans traits_10plus xys date
    clear ccData ccSignal conditionData experimentFolder
    
end


%% Part 4. plot bar graph of comparisons for each trait

cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
xSpacing = [0.86, 1.14; 1.86, 2.14; 2.86, 3.14];

% 1. tau
tau(1,1) = traits_steady{1,2}.means(6); % steady-state low
tau(1,2) = mean(traits_zero(:,6));      % only low signal (fluc)
tau(2,1) = traits_steady{1,3}.means(6); % steady-state average
tau(2,2) = mean(traits_50(:,6));        % 50% low signal (fluc)
tau(3,1) = traits_steady{1,4}.means(6); % steady-state high
tau(3,2) = mean(traits_100(:,6));       % only high signal (fluc)

tau_error(1,1) = traits_steady{1,2}.stds(6); % steady-state low
tau_error(1,2) = std(traits_zero(:,6));      % only low signal (fluc)
tau_error(2,1) = traits_steady{1,3}.stds(6); % steady-state average
tau_error(2,2) = std(traits_50(:,6));        % 50% low signal (fluc)
tau_error(3,1) = traits_steady{1,4}.stds(6); % steady-state high
tau_error(3,2) = std(traits_100(:,6));       % only high signal (fluc)

figure(1)
bar(tau*60)
hold on
errorbar(xSpacing,tau*60,tau_error*60,'.','Color',rgb('Black'))
legend('steady','fluc')
title('tau')
ylabel('tau')
plotName = strcat('I-60min-tau');
saveas(gcf,plotName,'epsc')
close(gcf)




% 2. division size
Vd(1,1) = traits_steady{1,2}.means(2); % steady-state low
Vd(1,2) = mean(traits_zero(:,2));      % only low signal (fluc)
Vd(2,1) = traits_steady{1,3}.means(2); % steady-state average
Vd(2,2) = mean(traits_50(:,2));        % 50% low signal (fluc)
Vd(3,1) = traits_steady{1,4}.means(2); % steady-state high
Vd(3,2) = mean(traits_100(:,2));       % only high signal (fluc)

Vd_error(1,1) = traits_steady{1,2}.stds(2); % steady-state low
Vd_error(1,2) = std(traits_zero(:,2));      % only low signal (fluc)
Vd_error(2,1) = traits_steady{1,3}.stds(2); % steady-state average
Vd_error(2,2) = std(traits_50(:,2));        % 50% low signal (fluc)
Vd_error(3,1) = traits_steady{1,4}.stds(2); % steady-state high
Vd_error(3,2) = std(traits_100(:,2));       % only high signal (fluc)

figure(2)
bar(Vd)
hold on
errorbar(xSpacing,Vd,Vd_error,'.','Color',rgb('Black'))
legend('steady','fluc')
title('Vd')
ylabel('division volume')
plotName = strcat('I-60min-Vd');
saveas(gcf,plotName,'epsc')
close(gcf)




% 3. mu
mu(1,1) = traits_steady{1,2}.means(7); % steady-state low
mu(1,2) = mean(traits_zero(:,7));      % only low signal (fluc)
mu(2,1) = traits_steady{1,3}.means(7); % steady-state average
mu(2,2) = mean(traits_50(:,7));        % 50% low signal (fluc)
mu(3,1) = traits_steady{1,4}.means(7); % steady-state high
mu(3,2) = mean(traits_100(:,7));       % only high signal (fluc)

mu_error(1,1) = traits_steady{1,2}.stds(7); % steady-state low
mu_error(1,2) = std(traits_zero(:,7));      % only low signal (fluc)
mu_error(2,1) = traits_steady{1,3}.stds(7); % steady-state average
mu_error(2,2) = std(traits_50(:,7));        % 50% low signal (fluc)
mu_error(3,1) = traits_steady{1,4}.stds(7); % steady-state high
mu_error(3,2) = std(traits_100(:,7));       % only high signal (fluc)

figure(3)
bar(mu)
hold on
errorbar(xSpacing,mu,mu_error,'.','Color',rgb('Black'))
legend('steady','fluc')
title('mu')
ylabel('mu')
plotName = strcat('I-60min-mu');
saveas(gcf,plotName,'epsc')
close(gcf)




% 4. added volume
added(1,1) = traits_steady{1,2}.means(3); % steady-state low
added(1,2) = mean(traits_zero(:,3));      % only low signal (fluc)
added(2,1) = traits_steady{1,3}.means(3); % steady-state average
added(2,2) = mean(traits_50(:,3));        % 50% low signal (fluc)
added(3,1) = traits_steady{1,4}.means(3); % steady-state high
added(3,2) = mean(traits_100(:,3));       % only high signal (fluc)

added_error(1,1) = traits_steady{1,2}.stds(3); % steady-state low
added_error(1,2) = std(traits_zero(:,3));      % only low signal (fluc)
added_error(2,1) = traits_steady{1,3}.stds(3); % steady-state average
added_error(2,2) = std(traits_50(:,3));        % 50% low signal (fluc)
added_error(3,1) = traits_steady{1,4}.stds(3); % steady-state high
added_error(3,2) = std(traits_100(:,3));       % only high signal (fluc)

figure(4)
bar(added)
hold on
errorbar(xSpacing,added,added_error,'.','Color',rgb('Black'))
legend('steady','fluc')
title('added volume')
ylabel('added volume')
plotName = strcat('I-60min-added');
saveas(gcf,plotName,'epsc')
close(gcf)




save('I-60min-data.mat','traits_steady','traits_zero','traits_50','traits_100')


%% Part 5. mann-whitney u test between bars of equal nutrient

clear
clc
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
load('I-60min-data.mat')

% loop through each trait (interdivision time, division volume, added volume)
trait_column = [6; 2; 3];

for tr = 1:length(trait_column)
    
    trait = trait_column(tr);
    
    % low
    low_steady = traits_steady{1,2}.values(:,trait); % steady-state low
    low_fluc = mean(traits_zero(:,trait));        % only low signal (fluc)
    [pL,hL] = ranksum(low_steady,low_fluc);
    [h_low,p_low] = ttest(low_steady,low_fluc);
    
    % ave
    ave_steady = traits_steady{1,3}.values(:,trait); % steady-state average
    ave_fluc = mean(traits_50(:,6));              % 50% low signal (fluc)
    [pA,hA] = ranksum(ave_steady,ave_fluc);
    [h_ave,p_ave] = ttest(ave_steady,ave_fluc);
    
    % high
    high_steady = traits_steady{1,4}.values(:,trait); % steady-state high
    high_fluc = mean(traits_100(:,trait));         % only high signal (fluc)
    [pH,hH] = ranksum(high_steady,high_fluc);
    [h_high,p_high] = ttest(high_steady,high_fluc);
    
    p_mann_whit(1,tr) = pL;
    p_mann_whit(2,tr) = pA;
    p_mann_whit(3,tr) = pH;
    
    p_studentsT(1,tr) = p_low;
    p_studentsT(2,tr) = p_ave;
    p_studentsT(3,tr) = p_high;

end


%% Part 6. plotting each nutrient group data over time

% goal: checking for stabilization
%       plot each trait for each over birth or div time

clear
clc
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
load('I-60min-data.mat')

%%
% loop through each trait (interdivision time, division volume, added volume)
trait_column = [6; 2; 3];

for tr = 1:length(trait_column)
    
    trait = trait_column(tr);
    
    % trait values
    low_trait = traits_zero(:,trait);
    ave_trait = traits_50(:,trait);
    high_trait = traits_100(:,trait);
    
    % times at birth
    low_births = traits_zero(:,4);
    ave_births = traits_50(:,4);
    high_births = traits_100(:,4);
    
    % times at birth
    low_divs = traits_zero(:,5);
    ave_divs = traits_50(:,5);
    high_divs = traits_100(:,5);
    
    figure(tr)
    scatter(low_births,low_trait,15,rgb('Indigo'),'filled')
    hold on
    scatter(ave_births,ave_trait,15,rgb('DarkGoldenRod'),'filled')
    hold on
    scatter(high_births,high_trait,15,rgb('DarkRed'),'filled')
    
    figure(tr+10)
    scatter(low_divs,low_trait,15,rgb('Indigo'),'filled')
    hold on
    scatter(ave_divs,ave_trait,15,rgb('DarkGoldenRod'),'filled')
    hold on
    scatter(high_divs,high_trait,15,rgb('DarkRed'),'filled')
    

end


