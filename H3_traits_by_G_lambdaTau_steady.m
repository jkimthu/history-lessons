%% Figure H3: traits by G=lambda*tau for steady environments

%  Goal: for a given experiment,
%
%        1.  calculate lambda*tau for all cell cycles after first 3h
%        2.  plot a scatter plot of trait values over lambda tau


%  Traits of interest: 
%
%       a) interdivision time
%       b) division size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score
%       f) added size




%  Last edit: jen, 2019 July 15
%  Commit: first commit, traits by lambda*tau for steady conditions

%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = 13:15; % use corresponding dataIndex values


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};


%% loop through steady conditions

for condition = 2:4  % low, ave and high only
    
    color = rgb(palette(condition));
    
    
    
    % Part 1. collect and concatenate final cell cycle data from each experiment
    
    % 1. for all experiments in dataset
    traits_all = [];
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
            %currentBinarySignal = binaryNutrientSignal(curveIDs == cellCyles_final(cc));
            %currentNscore = nScore(curveIDs == cellCyles_final(cc));
            currentTrack = tracks_final(cc);
            currentCC = cellCyles_final(cc);
            
            
            ccData(cc,1) = currentVolumes(1);       % V_birth
            ccData(cc,2) = currentVolumes(end);     % V_division
            ccData(cc,3) = currentVolumes(end) - currentVolumes(1); % added volume
            ccData(cc,4) = currentTimestamps(1);    % T_birth
            ccData(cc,5) = currentTimestamps(end);  % T_division
            ccData(cc,6) = currentTimestamps(end) - currentTimestamps(1); % interdivision time
            ccData(cc,7) = nanmean(currentGrowthRates); % mean growth rate
            ccData(cc,8) = nanstd(currentGrowthRates);  % standard deviation in growth rate
            %ccData(cc,9) = currentNscore(1);               % nScore
            
            %ccSignal{cc,1} = currentBinarySignal;        % binary signal (1 = high, 0 = low)
            ccData(cc,12) = currentTrack;
            ccData(cc,13) = currentCC;

            
            
        end
        clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly
        clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final
        clear binaryNutrientSignal nScore currentBinarySignal currentNscore
        clear isSwitch numShifts shiftStage currentTrack trackNum tracks_final currentCC
        clear bins fractions binned currentClass cl
        
        
        % 14. trim cell cycle data to avoid single point cell cycles
        addedVol = ccData(:,3);
        traits = ccData(addedVol > 0,:);
        %signal = ccSignal(addedVol > 0,:);
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
        %signal_final = signal(V_summed == 2,:);
        
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
    
    clear bubbletime ans traits_10plus xys date e experimentCount
    clear ccData ccSignal conditionData experimentFolder
    
    
    
    
    
    
    % Part 2. calculate lambda*tau and plot traits against it
    
    
    % 1. calculate lambda*tau for each cell cycle
    tau = traits_all(:,6);    % interdivision time (h)
    lambda = traits_all(:,7); % mean growth rate of cell cycle (1/h)
    G = lambda.*tau;
    
    
    
    % 2. isolate specific traits of interest
    tau = traits_all(:,6) * 60;
    Vd = traits_all(:,2);
    ratio = Vd./traits_all(:,1);
    mu = traits_all(:,7);
    added = traits_all(:,3);
    
    
    
    % 3. plot scatter with means binned by lambda*tau!
    sz = 7;
    
    figure(1)
    hold on
    scatter(G,tau,sz,color)
    
    figure(2)
    hold on
    scatter(G,Vd,sz,color)
    
    figure(3)
    hold on
    scatter(G,ratio,sz,color)
    
    figure(4)
    hold on
    scatter(G,mu,sz,color)
    
    figure(5)
    hold on
    scatter(G,added,sz,color)
    
    
end


% 5. label and save figures
dataFolder = '/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/H2_trait_vs_lambdaTau';
cd(dataFolder)

figure(1)
ylabel('interdivision time')
xlabel('lambda*tau')
xlim([0 1.8])
ylim([0 150])
title(strcat('H3-lambdaTau-tau-steady'))
plotName = strcat('H3-lambdaTau-tau-steady');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(2)
ylabel('division volume')
xlabel('lambda*tau')
xlim([0 1.8])
title(strcat('H3-lambdaTau-Vdiv-steady'))
plotName = strcat('H3-lambdaTau-Vdiv-steady');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(3)
ylabel('ratio')
xlabel('lambda*tau')
xlim([0 1.8])
title('H3-lambdaTau-ratio-steady')
plotName = strcat('H3-lambdaTau-ratio-steady');
saveas(gcf,plotName,'epsc')
close(gcf)

figure (4)
ylabel('mean growth rate (lambda)')
xlabel('lambda*tau')
xlim([0 1.8])
title('H3-lambdaTau-mu-steady')
plotName = strcat('H3-lambdaTau-mu-steady');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(5)
ylabel('added volume')
xlabel('lambda*tau')
xlim([0 1.8])
title('H3-lambdaTau-added-steady')
plotName = strcat('H3-lambdaTau-added-steady');
saveas(gcf,plotName,'epsc')
close(gcf)
    

    

