%% Figure H: traits by nClass -- continuously plotted

%  Goal: for a given experiment,
%
%        1.  classify all cell cycles after first 3h
%        2.  bin cell cycles by TYPE of nutrient classification
%               i.e. low to high, high to low, HLH, and LHL
%        3.  for each bin, plot a scatter plot of trait values over mean
%            nutrient score



%  Traits of interest: 
%
%       a) interdivision time
%       b) division size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score
%       f) added size




%  Last edit: jen, 2020 Jan 28
%  Commit: fix nScore and tau and measure size differences between signals
%          


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


% 0. initialize array for concatenation of final cell cycles from each experimental dataset
traits_all = [];


%% Part 1. collect and concatenate final cell cycle data from each experiment

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

clear bubbletime ans traits_10plus xys date
clear ccData ccSignal conditionData experimentFolder


%% Part 2. identify cell cycles by signal TYPE


% 1. bin traits by signal TYPE, based on classification
classifications = traits_all(:,14);
type = nan(length(classifications),1);

    % classes 2-5: high to low
    type(classifications == 2) = 1;
    type(classifications == 3) = 1;
    type(classifications == 4) = 1;
    type(classifications == 5) = 1;

    % classes 7-10: low to high
    type(classifications == 7) = 2;
    type(classifications == 8) = 2;
    type(classifications == 9) = 2;
    type(classifications == 10) = 2;

    % classes 11-14: high-low-high
    type(classifications == 11) = 3;
    type(classifications == 12) = 3;
    type(classifications == 13) = 3;
    type(classifications == 14) = 3;

    % classes 15-18: low-high-low
    type(classifications == 15) = 4;
    type(classifications == 16) = 4;
    type(classifications == 17) = 4;
    type(classifications == 18) = 4;

    % class 1 
    type(classifications == 1) = 5; % low only 
    
    % class 6
    type(classifications == 6) = 6;
    

%% Part 3. get sense of data by plotting Nscore and tau

tau = traits_all(:,6);
Nscore = traits_all(:,9);

% figure(1)
% scatter(Nscore,tau)
% xlabel('nScore'); ylabel('tau')
% 
% figure(2)
% histogram(Nscore,20)
% xlabel('Nscore'); ylabel('Count')
% 
% figure(3)
% histogram(tau,20)
% xlabel('Tau'); ylabel('Count')


% conclusion from figures 1-3:
%   most represented Nscore ~ 0.35-0.4      tau ~ 0.75-0.80
%                             0.30-0.35           0.70-0.75
%                             0.40-0.45           0.65-0.70

% make three plots, one for each of these three groups
% fix Nscore and tau, plot different traits across signal class
nScore_min = [0.35, 0.3, 0.4];
nScore_max = [0.399, 0.349, 0.449];

tau_min = [0.75, 0.70, 0.65];
tau_max = [0.799, 0.749, 0.699];

for gg = 1:length(tau_max)
    
    tau = traits_all(:,6);
    Nscore = traits_all(:,9);
    
    % isolate cell cycles with current nScore range
    current_n_min = nScore_min(gg);
    current_n_max = nScore_max(gg);
    
    traits_n1 = traits_all(Nscore > current_n_min,:);
    type_n1 = type(Nscore > current_n_min,:);
    Nscore_n1 = Nscore(Nscore > current_n_min);
    traits_n2 = traits_n1(Nscore_n1 <= current_n_max,:);
    type_n2 = type_n1(Nscore_n1 <= current_n_max,:);
    clear current_n_min current_n_max Nscore_n1 Nscore tau traits_n1
    clear type_n1
    
    
    % isolate cell cycles with current tau range
    tau_n = traits_n2(:,6);
    current_tau_min = tau_min(gg);
    current_tau_max = tau_max(gg);
    
    traits_t1 = traits_n2(tau_n > current_tau_min,:);
    type_t1 = type_n2(tau_n > current_tau_min,:);
    tau_n_t1 = tau_n(tau_n > current_tau_min);
    traits_t2 = traits_t1(tau_n_t1 <= current_tau_max,:);
    type_t2 = type_t1(tau_n_t1 <= current_tau_max,:);
    clear traits_t1 type_t1 traits_n2 tau_n_t1 tau_n
    clear current_tau_min current_tau_max
    
    % loop through signal types
    for tt = 1:6
        
        % 1. isolate traits from type of interest
        typeTraits = traits_t2(type_t2 == tt,:);
        
        % 2. isolate specific traits of interest
        Vd = typeTraits(:,2);
        added = typeTraits(:,3);
        ns = typeTraits(:,9);
        it = typeTraits(:,6);
        
        % 3. calculate mean and sem
        Vd_mean(tt) = nanmean(Vd);
        Vd_std(tt) = nanstd(Vd);
        Vd_count(tt) = length(Vd);
        Vd_sem(tt) = Vd_std(tt)./sqrt(Vd_count(tt));
        
        added_mean(tt) = nanmean(added);
        added_std(tt) = nanstd(added);
        added_count(tt) = length(added);
        added_sem(tt) = added_std(tt)./sqrt(added_count(tt));
        
        n_mean(tt) = nanmean(ns);
        n_std(tt) = nanstd(ns);
        n_sem(tt) = n_std(tt)./Vd_count(tt);
        
        it_mean(tt) =nanmean(it);
        it_std(tt) = nanstd(it);
        it_sem(tt) = it_std(tt)./Vd_count(tt);
        
    end
    
    % plot group data
%     figure(10+gg)
%     bar(Vd_mean)
%     hold on
%     er = errorbar([1,2,3,4,5,6],Vd_mean,Vd_sem,Vd_sem);
%     er.Color = [0 0 0];
%     er.LineStyle = 'none';
%     title(strcat('Group',num2str(gg)))
%     xlabel('nutrient signal')
%     ylabel('division v')
%     
%     figure(20+gg)
%     bar(added_mean)
%     hold on
%     er = errorbar([1,2,3,4,5,6],added_mean,added_sem,added_sem);
%     er.Color = [0 0 0];
%     er.LineStyle = 'none';
%     title(strcat('Group',num2str(gg)))
%     xlabel('nutrient signal')
%     ylabel('added v')
    
    figure(10+gg)
    bar(n_mean)
    hold on
    er = errorbar([1,2,3,4,5,6],n_mean,n_sem,n_sem);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title(strcat('Group',num2str(gg)))
    xlabel('nutrient signal')
    ylabel('nScore')
    ylim([0 0.5])
    
    figure(20+gg)
    bar(it_mean)
    hold on
    er = errorbar([1,2,3,4,5,6],it_mean,it_sem,it_sem);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title(strcat('Group',num2str(gg)))
    xlabel('nutrient signal')
    ylabel('tau')
    ylim([0 0.8])
   
    
end
%% Part 3. loop through signal TYPES, plotting continuous scatter for all traits

% type 1 = high to low
% type 2 = low to high
% type 3 = high-low-high
% type 4 = only high

for tt = 1:4
    
    

    

    
    
    % 3. bin traits by nScore for mean
    binSize = 0.5;
    numBins = 1/binSize;
    bins = ceil(nScore*numBins);
    
    binnedX = (1:max(bins))/numBins;
    binnedTau = accumarray(bins,tau,[],@mean);
    binnedVd = accumarray(bins,Vd,[],@mean);
    binnedRatio = accumarray(bins,ratio,[],@mean);
    binnedMu = accumarray(bins,mu,[],@mean);
    binnedAdded = accumarray(bins,added,[],@mean);

    
    
    % 4. plot scatter with means binned by nScore!
    figure(1)
    scatter(nScore,tau)
    hold on
    plot(binnedX,binnedTau)
    ylabel('interdivision time')
    xlabel('nScore')
    title(strcat('H-type-',num2str(t),'-tau'))
    plotName = strcat('H-type-',num2str(t),'-tau');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    scatter(nScore,Vd)
    hold on
    plot(binnedX,binnedVd)
    ylabel('division volume')
    xlabel('nScore')
    title(strcat('H-type-',num2str(t),'-Vdiv'))
    plotName = strcat('H-type-',num2str(t),'-Vdiv');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(3)
    scatter(nScore,ratio)
    hold on
    plot(binnedX,binnedRatio)
    ylabel('ratio')
    xlabel('nScore')
    title(strcat('H-type-',num2str(t),'-ratio'))
    plotName = strcat('H-type-',num2str(t),'-ratio');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(4)
    scatter(nScore,mu)
    hold on
    plot(binnedX,binnedMu)
    ylabel('growth rate')
    xlabel('nScore')
    title(strcat('H-type-',num2str(t),'-mu'))
    plotName = strcat('H-type-',num2str(t),'-mu');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(5)
    scatter(nScore,added)
    hold on
    plot(binnedX,binnedAdded)
    ylabel('added volume')
    xlabel('nScore')
    title(strcat('H-type-',num2str(t),'-added'))
    plotName = strcat('H-type-',num2str(t),'-added');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
end







