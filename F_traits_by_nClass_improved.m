%% Figure F: traits by nClass improved

%  Goal: for a given experiment,
%
%        1.  classify all cell cycles after first 3h
%        2.  bin cell cycles by nutrient classification
%        3.  plot mean and standard deviation in trait values across class



%  Traits of interest: 
%
%       a) interdivision time
%       b) division size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score




%  Last edit: jen, 2019 Mar 20
%  Commit: first commit, compiled 60 min data with grouped plotting


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
condition = 1;


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


% 0. initialize array for concatenation
traits_all = [];

%%

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
    
    
    
    % 5. isolate condition specific data
    color = rgb(palette(condition));
    
    
    
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
    
    
    
    % 9. calculate binary nutrient signals
    [binaryNutrientSignal, nScore] = nutrientScore(timescale,fullData);
    clear curveID growthRates growthRates_all
    
    
    
    % 10. isolate volume, isDrop, curveID, trackNum and timestamp data
    volumes = getGrowthParameter(fullData,'volume');
    isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
    curveIDs = getGrowthParameter(fullData,'curveFinder');
    trackNum = getGrowthParameter(fullData,'trackNum');
    timestamps_sec = getGrowthParameter(fullData,'timestamp');
    timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
    clear timestamps_sec
    
    
    % 11. identity unique cell cycles by ID number
    cellCycles = curveIDs(isDrop == 1);
    birthTimes = timestamps_hr(isDrop == 1);
    tracks = trackNum(isDrop == 1);
    clear fullData
    
    
    % 12. remove birth times prior to 3 hr
    birthTimes_post3 = birthTimes(birthTimes > 3);
    cellCycles_post3 = cellCycles(birthTimes > 3);
    tracks_post3 = tracks(birthTimes > 3);
    clear cellCycles birthTimes tracks
    
    
    % 13. remove birth times post bubbles
    cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
    tracks_final = tracks_post3(birthTimes_post3 < bubbletime(condition));
    clear cellCycles_post3 birthTimes_post3 tracks_post3
    
    
    % 14. for remaining cell cycles, identify:
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
    
    
    % 15. trim cell cycle data to avoid single point cell cycles
    addedVol = ccData(:,3);
    traits = ccData(addedVol > 0,:);
    signal = ccSignal(addedVol > 0,:);
    clear addedVol
    
    
    
    % 16. exclude outliers based on cell size (both birthsize and divsize)
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
    
    
    
    % 18. exclude interdivision times that are less than 10 min long
    % note: still need to justify exclusion of interdiv times < 10 visually
    interdivs = traits_final(:,6)*60;
    traits_10plus = traits_final(interdivs > 10,:);
    clear signal_final traits_final signal traits interdivs
    
    
    % 19. concatenate data across replicates
    traits_all = [traits_all; traits_10plus];
    size(traits_all)
    
end
    %%
    % some conditions do not have data after trimming
    
    % 19. bin traits by classification
    classifications = traits_all(:,14);
    
    tau = traits_all(:,6) * 60;
    Vd = traits_all(:,2);
    ratio = Vd./traits_all(:,1);
    mu = traits_all(:,7);
    nScore = traits_all(:,9);
    
    binnedCounts = accumarray(classifications,tau,[],@length);
    binnedTau_means = accumarray(classifications,tau,[],@mean);
    binnedTau_stds = accumarray(classifications,tau,[],@std);
    
    binnedVd_means = accumarray(classifications,Vd,[],@mean);
    binnedVd_stds = accumarray(classifications,Vd,[],@std);
    
    binnedRatio_means = accumarray(classifications,ratio,[],@mean);
    binnedRatio_stds = accumarray(classifications,ratio,[],@std);
    
    binnedMu_means = accumarray(classifications,mu,[],@mean);
    binnedMu_stds = accumarray(classifications,mu,[],@std);
    
    binnedNscore_means = accumarray(classifications,nScore,[],@mean);
    binnedNscore_stds = accumarray(classifications,nScore,[],@std);
    
    %%
    
    figure(1)
    bar(binnedTau_means,'FaceColor',color)
    hold on
    errorbar(binnedTau_means,binnedTau_stds,'.','Color',rgb('Black'))
    title('tau by signal class')
    
    figure(2)
    bar(binnedVd_means,'FaceColor',color)
    hold on
    errorbar(binnedVd_means,binnedVd_stds,'.','Color',rgb('Black'))
    title('division vol by signal class')
    
    figure(3)
    bar(binnedRatio_means,'FaceColor',color)
    hold on
    errorbar(binnedRatio_means,binnedRatio_stds,'.','Color',rgb('Black'))
    title('ratio by signal class')
    
    figure(4)
    bar(binnedMu_means,'FaceColor',color)
    hold on
    errorbar(binnedMu_means,binnedMu_stds,'.','Color',rgb('Black'))
    title('mu by signal class')
    
    figure(5)
    bar(binnedNscore_means,'FaceColor',color)
    hold on
    errorbar(binnedNscore_means,binnedNscore_stds,'.','Color',rgb('Black'))
    title('nScore by signal class')
    
    figure(6)
    bar(binnedCounts,'FaceColor',color)
    title('counts by signal class')
    
    
    % 19. save plots
    cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
    
    figure(1)
    plotName = strcat('E-fig1-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    figure(2)
    plotName = strcat('E-fig2-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    figure(3)
    plotName = strcat('E-fig3-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    figure(4)
    plotName = strcat('E-fig4-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    figure(5)
    plotName = strcat('E-fig5-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(6)
    plotName = strcat('E-fig6-all-c-',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    clear mu tau ratio Vb nScore
    
    %%
    
xClass = [1, 2,10, 3,9,11,17, 4,8,12,15, 5,7, 6];
xSpacing = [1, 3.8,4.8 7.6,8.6,9.6,10.6, 13.4,14.4,15.4,16.4, 19.2,20.2, 23];

    
figure(7)
xTau = binnedTau_means(xClass);
xTau_err = binnedTau_stds(xClass);
b = bar(xSpacing,xTau,0.8,'FaceColor',color);
hold on
errorbar(xSpacing,xTau,xTau_err,'.','Color',rgb('Black'))
ylabel('tau (min)')
title('interdivision times')
plotName = strcat('E-fig7-all-c-',num2str(condition));
saveas(gcf,plotName,'epsc')
close(gcf)


figure(8)
xVdiv = binnedVd_means(xClass);
xVdiv_err = binnedVd_stds(xClass);
b = bar(xSpacing,xVdiv,0.8,'FaceColor',color);
hold on
errorbar(xSpacing,xVdiv,xVdiv_err,'.','Color',rgb('Black'))
ylabel('Vdiv (um3)')
title('division volume')
plotName = strcat('E-fig8-all-c-',num2str(condition));
saveas(gcf,plotName,'epsc')
close(gcf)



figure(9)
xRatio = binnedRatio_means(xClass);
xRatio_err = binnedRatio_stds(xClass);
b = bar(xSpacing,xRatio,0.8,'FaceColor',color);
hold on
errorbar(xSpacing,xRatio,xRatio_err,'.','Color',rgb('Black'))
ylabel('Ratio')
title('Vdiv:Vbirth ratio')
plotName = strcat('E-fig9-all-c-',num2str(condition));
saveas(gcf,plotName,'epsc')
close(gcf)


figure(10)
xMu = binnedMu_means(xClass);
xMu_err = binnedMu_stds(xClass);
b = bar(xSpacing,xMu,0.8,'FaceColor',color);
hold on
errorbar(xSpacing,xMu,xMu_err,'.','Color',rgb('Black'))
ylabel('mu (1/h)')
title('growth rate (log2)')
plotName = strcat('E-fig10-all-c-',num2str(condition));
saveas(gcf,plotName,'epsc')
close(gcf)


figure(11)
xNscore = binnedNscore_means(xClass);
xNscore_err = binnedNscore_stds(xClass);
b = bar(xSpacing,xNscore,0.8,'FaceColor',color);
hold on
errorbar(xSpacing,xNscore,xNscore_err,'.','Color',rgb('Black'))
ylabel('Nscore')
title('average nutrient across cell cycle')
plotName = strcat('E-fig11-all-c-',num2str(condition));
saveas(gcf,plotName,'epsc')
close(gcf)





