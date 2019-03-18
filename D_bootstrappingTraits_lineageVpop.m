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




%  Last edit: jen, 2019 Mar 17
%  Commit: first commit, 


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
palette_below = {'LightSkyBlue','BlueViolet','Gold','LightCoral'};


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
        color = rgb(palette(condition));
        color_b = rgb(palette_below(condition));
        
        
        
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
        clear conditionData curveID growthRates growthRates_all
        
        
        
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
        divSize = traits_final(:,2);
        interDiv = traits_final(:,6);
        mu = traits_final(:,7);
        
        figure(1)
        qqplot(divSize)
        
        figure(2)
        qqplot(interDiv)
        
        figure(3)
        qqplot(mu)
        clear mu interDiv divSize
        
        
        % 18. bootstrap hypothesis testing!
        % note: still need to justify exclusion of interdiv times < 10 visually
        interdivs = traits_final(:,6)*60;
        traits_10plus = traits_final(interdivs > 10,:);
        signals_10plus = signal_final(interdivs > 10,:);
        clear signal_final traits_final signal traits
        
        classifications = traits_10plus(:,14);
        
        % some conditions do not have data after trimming
        if isempty(traits_10plus) == 1
            
            % filler
            
        else
            
            % A.  pull out a lineage with 4 or 5 cell cycles
            lineages = traits_10plus(:,12);
            uniqueLines = unique(lineages);
            uniqueCounts = hist(lineages,uniqueLines);
            
            % start with 5 consequtive, then 4
            longLines = uniqueLines(uniqueCounts == 5);
            numCC = 5;
            if isempty(longLines) == 1
                longLines = uniqueLines(uniqueCounts == 4);
                numCC = 4;
            end
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
                    disp(strcat('Cell cycles in lineage (',num2str(currenLine),') are not consecutive!'))
                end
                clear currentCycles isConsecutive
                
                
                % iii. generate random subsampling from entire population,
                %      repeat 100000 times
                repeats = 100000;
                subSampled_means = nan(5,repeats);
                
                for i = 1:repeats
                    
                    % collect subsampled population
                    sPop = nan(numCC,14);
                    for sample = 1:numCC
                        
                        sClass = currentClasses(sample);
                        sClass_idx = traits_10plus(classifications==sClass,:);
                        %sClass_signals = signal_10plus(classifications==sClass,:);
                        
                        row = randi(length(sClass_idx)); % random number generator from a uniform distribution of range 1 to length of sClass list
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
                    
                end
                
                save(strcat('D2_',date,'.mat'),'traits_10plus','signals_10plus','subSampled_means','exptData')
                
                % iv. determine probability of getting something more
                %     extreme than observed lineage (p-value)
                
                
                
            end
            

            
            
            % 18. plot PDF
            raw_vals = (1:max(assignedBins_raw)).*binSize; % raw division sizes

    
        end
        
        
        
    end
    
    % 19. save plots
    cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
    
    figure(1)
    plotName = strcat('C1-fig1-divSize-distributions-',date,'-binSize-',num2str(binSize));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    figure(2)
    plotName = strcat('C1-fig2-traitComparisons-1stdAboveorBelow-',date,'-c1');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(3)
    plotName = strcat('C1-fig3-traitComparisons-1stdAboveorBelow-',date,'-c2');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(4)
    plotName = strcat('C1-fig4-traitComparisons-1stdAboveorBelow-',date,'-c3');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(5)
    plotName = strcat('C1-fig5-traitComparisons-1stdAboveorBelow-',date,'-c4');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    clc
    
    
end


