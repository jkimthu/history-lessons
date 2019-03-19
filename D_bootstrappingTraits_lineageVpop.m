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




%  Last edit: jen, 2019 Mar 19
%  Commit: test statistic is standard deviation of 5 traits, lineage length
%  is 4 cell cycles


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
condition = 1;

% 0. initialize number of loops for bootstrapping
repeats = 10000;


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
    ];

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
    

    
    % 3. load measured experiment data    
    %experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    %cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    % 4. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType filename
    
    
    
    
    % 5. isolate condition specific data
    conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
    
    
    
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
    clear conditionData curveID growthRates growthRates_all timescale
    
    
    
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
        bins_percent = ceil(fractions'*5);
        binned = round(accumarray(bins_percent,currentBinarySignal,[],@mean));
        
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
    clear ccData addedVol ccSignal
    
    
    
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
    
    
    
    % 17. qq plots to determine whether data has normal distribution
%     divSize = traits_final(:,2);
%     interDiv = traits_final(:,6);
%     mu = traits_final(:,7);
%     
%     figure(1)
%     qqplot(divSize)
%     title('QQ plot of division size')
%     plotName = strcat('D-QQ-divSize-',date,'-c',num2str(condition));
%     saveas(gcf,plotName,'epsc')
%     close(gcf)
%         
%     figure(2)
%     qqplot(interDiv)
%     title('QQ plot of interdivision time')
%     plotName = strcat('D-QQ-tau-',date,'-c',num2str(condition));
%     saveas(gcf,plotName,'epsc')
%     close(gcf)
%     
%     figure(3)
%     qqplot(mu)
%     title('QQ plot of growth rate')
%     plotName = strcat('D-QQ-mu-',date,'-c',num2str(condition));
%     saveas(gcf,plotName,'epsc')
%     close(gcf)
%     clear mu interDiv divSize
    
    
    % 18. bootstrap hypothesis testing!
    % note: still need to justify exclusion of interdiv times < 10 visually
    interdivs = traits_final(:,6)*60;
    traits_10plus = traits_final(interdivs > 10,:);
    signals_10plus = signal_final(interdivs > 10,:);
    clear signal_final traits_final signal traits interdivs
    
    % plot distribution of nutrient signal classifications
    classifications = traits_10plus(:,14);
    figure(1)
    hist(classifications,15)
    title(strcat('Histogram of nutrient signal classifications :',date))
    ylabel('Frequency')
    xlabel('Signal class')
    xlim([0 18])
    plotName = strcat('D-histogram-nClass-',date,'-c',num2str(condition));
    saveas(gcf,plotName,'epsc')
    close(gcf)

    
    % some conditions do not have data after trimming
    if isempty(traits_10plus) == 1
        
        continue
        
    else
        
        % A.  pull out a lineage with 4 or 5 cell cycles
        lineages = traits_10plus(:,12);
        uniqueLines = unique(lineages);
        uniqueCounts = hist(lineages,uniqueLines);
        
        % start with 5 consequtive, then 4
        numCC = 4;
        longLines = uniqueLines(uniqueCounts == numCC);
        clear uniqueLines uniqueCounts
        
        
        % B.  use bootstrapping hypothesis testing to determine
        %     likelihood for each long lineage
        for l = 1:length(longLines)
            
            
            % i. determine traits for each cycle within lineage
            currentLine = longLines(l);
            currentTraits = traits_10plus(lineages == currentLine,:);
            currentClasses = currentTraits(:,14);
            
            
            % ii. if lineage cell cycles are not consecutive, note this!
            currentCycles = currentTraits(:,13);
            isConsecutive = diff(currentCycles);
            if mean(isConsecutive) ~= 1
                disp(strcat('Cell cycles in lineage (',num2str(currentLine),') are not consecutive!'))
            end
            clear currentCycles isConsecutive
        
            
            % iii. generate random subsampling from entire population,
            %      repeat 10000 times
            subSampled_means = nan(repeats,5);
            subSampled_stds = nan(repeats,5);
            for i = 1:repeats
                
                % collect subsampled population
                sPop = nan(numCC,14);
                for sample = 1:numCC
                    
                    sClass = currentClasses(sample);
                    sClass_idx = traits_10plus(classifications==sClass,:);
                    [sClass_list,~] = size(sClass_idx);
                    
                    row = randi(sClass_list); % random number generator from a uniform distribution of range 1 to length of sClass list
                    sPop(sample,:) =  sClass_idx(row,:);
                    
                end
                clear row sClass_idx sClass sample
                
                % average test statistics
                sPop_mean = mean(sPop);
                subSampled_means(i,1) = sPop_mean(6)*60; % interdivision time (min)
                subSampled_means(i,2) = sPop_mean(1); % birth size
                subSampled_means(i,3) = sPop_mean(2)/sPop_mean(1); % V_div/V_birth ratio
                subSampled_means(i,4) = sPop_mean(7); % mean growth rate of cell cycle
                subSampled_means(i,5) = sPop_mean(9); % nutrient score
                
                sPop_std = std(sPop);
                subSampled_stds(i,1) = sPop_std(6)*60;
                subSampled_stds(i,2) = sPop_std(1); % birth size
                subSampled_stds(i,3) = sPop_std(2)/sPop_std(1); % V_div/V_birth ratio
                subSampled_stds(i,4) = sPop_std(7); % mean growth rate of cell cycle
                subSampled_stds(i,5) = sPop_std(9); % nutrient score
                
            end
            clear i sPop sPop_mean
            
            
            % iv. histograms and means of bootstrapped test stats
            for h = 1:5
                
                meanStat(1,h) = mean(subSampled_means(:,h));
                stdStat(1,h) = mean(subSampled_stds(:,h));
                
                figure(1)
                subplot(1,5,h)
                hist(subSampled_means(:,h))
                title(num2str(meanStat(1,h)))
                
                figure(2)
                subplot(1,5,h)
                hist(subSampled_stds(:,h))
                title(num2str(stdStat(1,h)))
                
            end
            figure(1)
            plotName = strcat('D-subsample-histograms-mean-',date,'-c',num2str(condition),'-line',num2str(l));
            saveas(gcf,plotName,'epsc')
            close(gcf)
            clear h plotName
            
            figure(2)
            plotName = strcat('D-subsample-histograms-std-',date,'-c',num2str(condition),'-line',num2str(l));
            saveas(gcf,plotName,'epsc')
            close(gcf)
            clear h plotName
            
            
            % v. determine probability of getting something more
            %    extreme than observed lineage (p-value)
            lineage_means(l,1) = mean(currentTraits(:,6)) * 60;
            lineage_means(l,2) = mean(currentTraits(:,1));
            lineage_means(l,3) = mean(currentTraits(:,2)./currentTraits(:,1));
            lineage_means(l,4) = mean(currentTraits(:,7));
            lineage_means(l,5) = mean(currentTraits(:,9));
            
            lineage_stds(l,1) = std(currentTraits(:,6)) * 60;
            lineage_stds(l,2) = std(currentTraits(:,1));
            lineage_stds(l,3) = std(currentTraits(:,2)./currentTraits(:,1));
            lineage_stds(l,4) = std(currentTraits(:,7));
            lineage_stds(l,5) = std(currentTraits(:,9));
            
            % more extreme greater than distance between lineage mean
            % and bootstrapped mean
            distance = abs(stdStat - lineage_stds(l,:));
            
            
            % calculate p-vals for each test statistic
            extreme_lows = stdStat - distance;
            extreme_highs = stdStat + distance;
            
            for col = 1:5
                count_lows(:,col) = subSampled_stds(:,col) < extreme_lows(col);
                count_highs(:,col) = subSampled_stds(:,col) > extreme_highs(col);
            end
            
            pVals(l,:) = (sum(count_lows)+sum(count_highs))./repeats;
            
        end
        clear extreme_lows extreme_highs count_lows count_highs distance l
        clear l longLines subSampled_means meanStat
        
    end
    
    save(strcat('D-',date,'-std-c1-length4'),'pVals','lineage_stds')
    
    clear signals_10plus traits_10plus classifications
    clear pVals lineage_means interdivs
    clear currentClasses currentLine currentTraits numCC
end
clear specificColumn bubbletime repeats

%%

directory = dir(strcat('D*.mat'));
names = {directory.name};
trait_name = {'tau','Vbirth','ratio','mu','nScore'};


for i = 1:length(names)
    
    load(names{i})

    figure()
    for trait = 1:5
        
        subplot(1,5,trait)
        hist(pVals(:,trait))
        title(trait_name{trait})
        xlim([0 1])
        
        if trait == 3
            xlabel(names{i})
        end
        
    end
    plotName = strcat('D-pVals-directory-',num2str(i));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
end
