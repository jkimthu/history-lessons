%% Figure D: bootstrapping to determine whether lineage history influences cell cycle traits


%  Goal: for a given experiment,
%
%        1.  pull out a lineage with maximum number of cell cycles
%        2.  determine traits (including nutrient history) for each cycle
%        3.  use bootstrapping hypothesis testing to re-sample population
%        4.  determine probability of getting something more extreme than
%            observed lineage (p-value)
%        5.  population subsample should take nutrient scores of lineage
%            sample into account



%  Traits of interest: 
%
%       a) interdivision time
%       b) birth size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score
%       f) added volume



%  Last edit: jen, 2019 July 24
%  Commit: first commit, bootstrap test on steady average to determine
%          lineage dependencies


%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 0. define method of calculating growth rate
specificColumn = 3;      % for selecting log2 column in growthRates


% 0. initialize experiments and condition to use in analysis
exptArray = 13:15; % use corresponding dataIndex values
condition = 3;

% 0. initialize number of generations in lineage
numCC = 4;

% 0. initialize number of loops for bootstrapping
repeats = 10000;


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
    %experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    %cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    % 4. compile experiment data matrix
    xy_start = xys(condition,1);
    xy_end = xys(condition,end);
    conditionData = buildDM(D5,T,xy_start,xy_end,index,expType);
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
    
    
    
    % 8. isolate volume, isDrop, curveID, trackNum and timestamp data
    volumes = getGrowthParameter(fullData,'volume');
    isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
    curveIDs = getGrowthParameter(fullData,'curveFinder');
    trackNum = getGrowthParameter(fullData,'trackNum');
    timestamps_sec = getGrowthParameter(fullData,'timestamp');
    timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
    clear timestamps_sec
    
    
    % 9. identity unique cell cycles by ID number
    cellCycles = curveIDs(isDrop == 1);
    birthTimes = timestamps_hr(isDrop == 1);
    tracks = trackNum(isDrop == 1);
    clear fullData
    
    
    % 10. remove birth times prior to 3 hr
    birthTimes_post3 = birthTimes(birthTimes > 3);
    cellCycles_post3 = cellCycles(birthTimes > 3);
    tracks_post3 = tracks(birthTimes > 3);
    clear cellCycles birthTimes tracks
    
    
    % 11. remove birth times post bubbles
    cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
    tracks_final = tracks_post3(birthTimes_post3 < bubbletime(condition));
    clear cellCycles_post3 birthTimes_post3 tracks_post3
    
    
    % 12. for remaining cell cycles, identify:
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
        
        ccData(cc,12) = currentTrack;
        ccData(cc,13) = currentCC;
        
    end
    clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly
    clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final
    clear binaryNutrientSignal nScore currentBinarySignal currentNscore
    clear isSwitch numShifts shiftStage currentTrack trackNum tracks_final currentCC
    clear bins fractions binned currentClass cl
    
    
    % 13. trim cell cycle data to avoid single point cell cycles
    addedVol = ccData(:,3);
    traits = ccData(addedVol > 0,:);
    clear ccData addedVol ccSignal
    
    
    
    % 14. exclude outliers based on cell size (both birthsize and divsize)
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
    clear switches
    
    
    
    % 15. bootstrap hypothesis testing!
    % note: still need to justify exclusion of interdiv times < 10 visually
    interdivs = traits_final(:,6)*60;
    traits_10plus = traits_final(interdivs > 10,:);
    clear signal_final traits_final signal traits interdivs
    
    
    % some conditions do not have data after trimming
    if isempty(traits_10plus) == 1
        
        continue
        
    else
        
        % A.  pull out a lineage with 4 or 5 cell cycles
        lineages = traits_10plus(:,12);
        uniqueLines = unique(lineages);
        uniqueCounts = hist(lineages,uniqueLines);

        longLines = uniqueLines(uniqueCounts == numCC);
        clear uniqueLines uniqueCounts
        
        
        % B.  use bootstrapping hypothesis testing to determine
        %     likelihood for each long lineage
        for ll = 1:length(longLines)
            
            
            % i. determine traits for each cycle within lineage
            currentLine = longLines(ll);
            currentTraits = traits_10plus(lineages == currentLine,:);

            
            % ii. if lineage cell cycles are not consecutive, note this!
            currentCycles = currentTraits(:,13);
            isConsecutive = diff(currentCycles);
            if mean(isConsecutive) ~= 1
                disp(strcat('Cell cycles in lineage (',num2str(currentLine),') are not consecutive!'))
            end
            clear currentCycles isConsecutive currentLine
        
            
            % iii. generate random subsampling from entire population,
            %      repeat 10000 times
            subSampled_means = nan(repeats,5);
            subSampled_stds = nan(repeats,5);
            for ii = 1:repeats
                
                % collect subsampled population
                sPop = nan(numCC,13);
                for sample = 1:numCC
                    
                    % in steady, all cell cycles are fair game
                    % sample with replacement
                    cc_total = length(traits_10plus);
                    
                    row = randi(cc_total); % random number generator from a uniform distribution of range 1 to length of sClass list
                    sPop(sample,:) =  traits_10plus(row,:);
                    
                end
                clear row cc_total sClass sample
                
                % average test statistics
                sPop_mean = mean(sPop);
                subSampled_means(ii,1) = sPop_mean(6)*60; % interdivision time (min)
                subSampled_means(ii,2) = sPop_mean(1); % birth size
                subSampled_means(ii,3) = sPop_mean(2)/sPop_mean(1); % V_div/V_birth ratio
                subSampled_means(ii,4) = sPop_mean(7); % mean growth rate of cell cycle
                subSampled_means(ii,5) = sPop_mean(9); % nutrient score
                subSampled_means(ii,6) = sPop_mean(3); % added vol
                
                sPop_std = std(sPop);
                subSampled_stds(ii,1) = sPop_std(6)*60;
                subSampled_stds(ii,2) = sPop_std(1); % birth size
                subSampled_stds(ii,3) = sPop_std(2)/sPop_std(1); % V_div/V_birth ratio
                subSampled_stds(ii,4) = sPop_std(7); % mean growth rate of cell cycle
                subSampled_stds(ii,5) = sPop_std(9); % nutrient score
                subSampled_stds(ii,6) = sPop_std(3); % added vol
                
            end
            clear i sPop sPop_mean
            
            
            % iv. histograms and means of bootstrapped test stats
            for h = 1:6
                
                meanStat(1,h) = mean(subSampled_means(:,h));
                stdStat(1,h) = mean(subSampled_stds(:,h));
                
            end
            
            
            % v. determine probability of getting something more
            %    extreme than observed lineage (p-value)
            lineage_means(ll,1) = mean(currentTraits(:,6)) * 60;
            lineage_means(ll,2) = mean(currentTraits(:,1));
            lineage_means(ll,3) = mean(currentTraits(:,2)./currentTraits(:,1));
            lineage_means(ll,4) = mean(currentTraits(:,7));
            lineage_means(ll,5) = mean(currentTraits(:,9));
            lineage_means(ll,6) = mean(currentTraits(:,3));
            
            lineage_stds(ll,1) = std(currentTraits(:,6)) * 60;
            lineage_stds(ll,2) = std(currentTraits(:,1));
            lineage_stds(ll,3) = std(currentTraits(:,2)./currentTraits(:,1));
            lineage_stds(ll,4) = std(currentTraits(:,7));
            lineage_stds(ll,5) = std(currentTraits(:,9));
            lineage_stds(ll,6) = std(currentTraits(:,3));
            
            % more extreme greater than distance between lineage mean
            % and bootstrapped mean
            distance = abs(stdStat - lineage_stds(ll,:));
            
            
            % calculate p-vals for each test statistic
            extreme_lows = stdStat - distance;
            extreme_highs = stdStat + distance;
            
            for col = 1:5
                count_lows(:,col) = subSampled_stds(:,col) < extreme_lows(col);
                count_highs(:,col) = subSampled_stds(:,col) > extreme_highs(col);
            end
            
            pVals(ll,:) = (sum(count_lows)+sum(count_highs))./repeats;
            
        end
        clear extreme_lows extreme_highs count_lows count_highs distance l
        clear l longLines subSampled_means meanStat
        
    end
    
    save(strcat('D-',date,'-std-c',num2str(condition),'-length',num2str(numCC)),'pVals','lineage_stds','lineage_means')
    
    clear signals_10plus traits_10plus classifications
    clear pVals lineage_means interdivs lineage_stds lineage_stds
    clear currentClasses currentLine currentTraits
end
clear specificColumn bubbletime repeats

%%

directory = dir(strcat('D*.mat'));
names = {directory.name};
trait_name = {'tau','Vbirth','ratio','mu','nScore'};


for ii = 1:length(names)
    
    load(names{ii})

    figure()
    for trait = 1:5
        
        subplot(1,5,trait)
        hist(pVals(:,trait))
        title(trait_name{trait})
        xlim([0 1])
        
        if trait == 3
            xlabel(names{ii})
        end
        
    end
    plotName = strcat('D-pVals-directory-',num2str(ii));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
end
