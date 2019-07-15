%% Figure H: traits by G=lambda*tau 

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
%  Commit: edit to control color scheme, don't bin lambda*tau

%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
%experimentCount = length(dataIndex);


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = 13:15; % use corresponding dataIndex values
condition = 1;


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

clear bubbletime ans traits_10plus xys date e experimentCount
clear ccData ccSignal conditionData experimentFolder


%% Part 2. calculate and identify cell cycles by lambda*tau

clear dataIndex specificColumn specificGrowthRate

% 1. calculate lambda*tau for each cell cycle
tau = traits_all(:,6);    % interdivision time (h)
lambda = traits_all(:,7); % mean growth rate of cell cycle (1/h)
G = lambda.*tau;


% 2. identify cell cycles by course signal classifications
classifications = traits_all(:,14);
type = nan(length(classifications),1);

    % classes 1 and 6 are full high and full low
    type(classifications == 1) = 5; % only low
    type(classifications == 6) = 6; % only high
    
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

    
    
% 3. confirm all shifting signals are accounted for
noShifts = find(isnan(type));

if isempty(noShifts) == 0
    error('Not all shifting signals are accounted for!')
else
    disp('All shifting signals are accounted for!')
end


%% Part 3. loop through signal TYPES, plotting continuous scatter for all traits


clear allHigh allLow classifications lambda tau classRules index exptArray condition


% 0. initialize colors for plotting
palette = {'SteelBlue','Chocolate','DarkOrchid','MediumSeaGreen','Black','SlateGray'};

for t = 1:6
    
    % type 1 = high to low
    % type 2 = low to high
    % type 3 = high-low-high
    % type 4 = low-high-low
    % type 5 = only low
    % type 6 = only high
    
    
    
    % 1. isolate traits from type of interest
    typeTraits = traits_all(type == t,:);
    type_Gs = G(type == t);
    
    
    % 2. isolate specific traits of interest
    tau = typeTraits(:,6) * 60;
    Vd = typeTraits(:,2);
    ratio = Vd./typeTraits(:,1);
    mu = typeTraits(:,7);
    added = typeTraits(:,3);
    
    
    
    % 3. plot scatter with means binned by lambda*tau!
    color = rgb(palette(t));
    sz = 7;
    
    figure(1)
    hold on
    scatter(type_Gs,tau,sz,color)

    figure(2)
    hold on
    scatter(type_Gs,Vd,sz,color)
    
    figure(3)
    hold on
    scatter(type_Gs,ratio,sz,color)
  
    figure(4)
    hold on
    scatter(type_Gs,mu,sz,color)

    figure(5)
    hold on
    scatter(type_Gs,added,sz,color)
    
end


% 5. label and save figures
dataFolder = '/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/H2_trait_vs_lambdaTau';
cd(dataFolder)

for ff = 1:5
    
    figure(1)
    ylabel('interdivision time')
    xlabel('lambda*tau')
    xlim([0 1.8])
    title(strcat('H2-lambdaTau-tau'))
    plotName = strcat('H2-lambdaTau-tau');
    saveas(gcf,plotName,'epsc')
    %close(gcf)
    
    
    figure(2)
    ylabel('division volume')
    xlabel('lambda*tau')
    xlim([0 1.8])
    title(strcat('H2-lambdaTau-Vdiv'))
    plotName = strcat('H2-lambdaTau-Vdiv');
    saveas(gcf,plotName,'epsc')
    %close(gcf)
    
    figure(3)
    ylabel('ratio')
    xlabel('lambda*tau')
    xlim([0 1.8])
    title(strcat('H2-lambdaTau-ratio'))
    plotName = strcat('H2-lambdaTau-ratio');
    saveas(gcf,plotName,'epsc')
    %close(gcf)
    
    figure (4)
    ylabel('mean growth rate (lambda)')
    xlabel('lambda*tau')
    xlim([0 1.8])
    title(strcat('H2-lambdaTau-mu'))
    plotName = strcat('H2-lambdaTau-mu');
    saveas(gcf,plotName,'epsc')
    %close(gcf)
    
    figure(5)
    ylabel('added volume')
    xlabel('lambda*tau')
    xlim([0 1.8])
    title(strcat('H2-lambdaTau-added'))
    plotName = strcat('H2-lambdaTau-added');
    saveas(gcf,plotName,'epsc')
    %close(gcf)
end



