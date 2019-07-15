%% Figure D2: autocorrelation analysis of traits between generations


%  Goal: for a given experiment,
%
%        1.  pull out a lineage with a set number of cell cycles
%        2.  determine traits (including nutrient history) for each cycle
%        3.  vary the autocorrelation shift
%        4.  comupte autocorrelation across shifts possible
%        5.  color points based on signal classifications of two
%            generations tested



%  Traits of interest: 
%
%       a) interdivision time
%       b) division size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score
%       f) added size
%       g) birth size




%  Last edit: jen, 2019 July 15
%  Commit: correlations between mom-daughter and grandma-granddaughter pairs
%          across all timescales, store values and plot


%  OK let's go!

%% data processing and autocorrelation analysis

% PART 0. initialize analysis

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




% 0. initalize all possible condition types to loop through
ts_all = {2:4;5:7;9:12;13:15};

for ts = 1:length(ts_all)
    

    % 0. initialize array of experiments to use in analysis, then loop through each
    exptArray = ts_all{ts}; % use corresponding dataIndex values
    
  
    
    % PART 1. collect and concatenate cell cycle data from each fluctuating experiment
    
    for condition = 1:4 % condition goes first to accumulate all same condition data across replicates
        
        
        % 0. initialize array for concatenation of final cell cycles from each experimental dataset
        traits_all = [];
        
        
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
                
                
%                 isSwitch = [0; diff(currentBinarySignal)] ~= 0;
%                 numShifts = sum(isSwitch);
%                 if numShifts > 0
%                     switches{cc,1} = find(isSwitch == 1);
%                     shiftStage = switches{cc}(1)/length(currentBinarySignal); % fraction into cell cycle of first switch
%                 else
%                     switches{cc,1} = NaN;
%                     shiftStage = NaN;
%                 end
                
                
                %ccData(cc,10) = numShifts;
                %ccData(cc,11) = shiftStage; % first if multiple
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
            
            clear birthSize_median birthSize_std divSize_median divSize_std
            clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
            clear V_birth_binary V_division_binary birth_outliers div_outliers V_div V_birth V_summed interdiv
            clear switches
            
            
            
            % 16. exclude interdivision times that are less than 10 min long
            % note: still need to justify exclusion of interdiv times < 10 visually
            interdivs = traits_final(:,6)*60;
            traits_10plus = traits_final(interdivs > 10,:);
            clear traits_final  traits interdivs
            
            
            % 17. concatenate data across replicates
            traits_all = [traits_all; traits_10plus];
            size(traits_all)
            
        end
        
        ccData_compiled{ts,condition} = ccData;
        clear bubbletime ans traits_10plus xys date
        clear ccData conditionData experimentFolder
        clear index e
        
        
        
        
        
        % PART 2. mom-daughter autocorrelation
        
        
        % 0. determine minimum cell cycles per lineage
        minCC = 2;
        
        
        % 1. pull out lineages of at least min cell cycles
        lineages = traits_all(:,12);
        uniqueLines = unique(lineages);
        uniqueCounts = hist(lineages,uniqueLines);
        longLines = uniqueLines(uniqueCounts >= minCC);
        
        
        % 2. generate n (mom) and n+1 (daughter) trait matrices
        mom_traits = [];
        daughter_traits = [];
        
        for l = 1:length(longLines)
            
            counter = 1;
            
            % i. isolate traits for each cycle within lineage
            currentLine = longLines(l);
            currentTraits = traits_all(lineages == currentLine,:);
            currentCCs = sum(lineages == currentLine);
            
            
            % ii . based on number of cell cycles in current line, assign
            while counter < currentCCs
                
                mom_traits = [mom_traits; currentTraits(counter,:)];
                daughter_traits = [daughter_traits; currentTraits(counter+1,:)];
                counter = counter + 1;
                
            end
            
            
        end
        clear l counter currentCCs currentTraits currentLine
        
        
        
        % 3. plot n+1 (daughter) vs n (mom) for traits of interest
        toi = {'tau','Vd','Vb','mu','added','nScore'};
        toi_col = [6; 2; 1; 7; 3; 9];
        
        for t = 1:length(toi_col)
            
            
            % i. isolate trait data
            t_mom = mom_traits(:,toi_col(t));
            t_daughter = daughter_traits(:,toi_col(t));
            
            
            % ii. fit linear regression
            p = polyfit(t_mom,t_daughter,1);
            x = t_mom;
            y = p(1)*x + p(2);
            
            
            % iii. calculate correlation coefficient
            r = corrcoef(t_mom,t_daughter);
            R1(t) = r(1,2);
            clear p x y r txt t_daughter t_mom
            
        end
        
        R1_compiled{condition}(ts,:) = R1;
        clear mom_traits daughter_traits
        clear t
        
        
        
        
        % PART 3. grandma-granddaughter autocorrelation
        
        % 0. determine minimum cell cycles per lineage
        minCC_shift2 = 3;
        
        
        % 1. pull out lineages of at least min cell cycles
        longerLines = uniqueLines(uniqueCounts >= minCC_shift2);
        
        
        % 2. generate n (mom) and n+1 (daughter) trait matrices
        gma_traits = [];
        gda_traits = [];
        
        for l = 1:length(longerLines)
            
            counter2 = 1;
            
            % i. isolate traits for each cycle within lineage
            currentLine = longerLines(l);
            currentTraits = traits_all(lineages == currentLine,:);
            currentCCs = sum(lineages == currentLine);
            
            % ii . based on number of cell cycles in current line, assign
            while counter2 < (currentCCs-1)
                
                gma_traits = [gma_traits; currentTraits(counter2,:)];
                gda_traits = [gda_traits; currentTraits(counter2+2,:)];
                counter2 = counter2 + 1;
                
            end
            
        end
        clear l counter currentCCs currentTraits currentLine
        
        
        
        % 3. plot n+1 (daughter) vs n (mom) for traits of interest
        for t = 1:length(toi_col)
            
            
            % i. isolate trait data
            t_granny = gma_traits(:,toi_col(t));
            t_grandkid = gda_traits(:,toi_col(t));
            
            
            % ii. fit linear regression
            p = polyfit(t_granny,t_grandkid,1);
            x = t_granny;
            y = p(1)*x + p(2);
            
            
            % iii. calculate correlation coefficient
            r = corrcoef(t_granny,t_grandkid);
            R2(t) = r(1,2);
            clear p x y r txt t_granny t_grandkid
            
        end
        
        R2_compiled{condition}(ts,:) = R2;
        clear mom_traits daughter_traits
        clear t
        
    end
end

saveFolder = '/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/D2_autocorrelations';
cd(saveFolder)

save('autocorrelations_compiled','R1_compiled','R2_compiled','ccData_compiled','toi')

%% data visualizaton

% goal: plot autocorrelation function,
% where x axis is the number of generations and y axis is correlation coefficient (R)


% 0. initialize analysis
clear
clc 
dataFolder = '/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/D2_autocorrelations';
cd(dataFolder)

load('autocorrelations_compiled.mat')



% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};



% 1. for each trait of interest, collect data from all timescales
for tr = 1:length(toi)
    
    % 2. loop through each condition
    for cond = 1:4
        
        currentR1 = R1_compiled{cond}(:,tr);
        currentR2 = R2_compiled{cond}(:,tr);
        currentR0 = [1;1;1;1];
        
        
        % 3. concatenate data such that column = generations apart; row = timescale
        compiledR = [currentR0, currentR1, currentR2];
        generation = [0,1,2];
        
        
        % 4. determine plotting color and shade
        color = rgb(palette(cond));
        lineType = {'-', '--', '-.', ':'};
        
        
        % 5. plot autocorrelation function for each timescale
        for ts = 1:4
           
            figure(tr)
            hold on
            plot(generation,compiledR(ts,:),'Color',color,'LineStyle',lineType{ts})
            
            if ts == 4
                title(strcat('autocorrelation -',toi{tr}))
                xlabel('generations apart')
                ylabel('correlation coefficient (R)')
            end
            
        end
        
    end
    
    % 6. save plot in active folder
    plotName = strcat('D2-autocorr-',toi{tr});
    saveas(gcf,plotName,'epsc')
    
    %close(gcf)
    %clc
    
end

