%% Figure C1: traits associated with being above or below average _ division size

%  Goal: for a given experiment,
%
%        1.  determine the mean size at division
%        2.  plot the distribution
%        3.  split condition data into two populations: above and below average
%        4.  compare the mean value for a variety of traits


%  Traits of interest: 
%
%       a) interdivision time
%       b) birth size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle / mean growth rate of population
%       e) nutrient score
%       f) cell cycle fraction at the point of a shift (shift stage)



%  Last edit: jen, 2019 Mar 8
%  Commit: first commit, not yet including nutrient score or shift stage


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
exptArray = 10:15; % use corresponding dataIndex values


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
palette_below = {'LightSkyBlue','BlueViolet','Gold','LightCoral'};



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
        
                
        %nData = nutrientScore(timescale,fullData);
        [binaryNutrientSignal, nScore] = nutrientScore(timescale,fullData);
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
        ccData = nan(length(cellCyles_final),11);
        
        for cc = 1:length(cellCyles_final)
            
            % isolate data for current cell cycle
            currentVolumes = volumes(curveIDs == cellCyles_final(cc));
            currentTimestamps = timestamps_hr(curveIDs == cellCyles_final(cc));
            currentGrowthRates = growthRates_fullOnly(curveIDs == cellCyles_final(cc));
            currentBinarySignal = binaryNutrientSignal(curveIDs == cellCyles_final(cc));
            currentNscore = nScore(curveIDs == cellCyles_final(cc));
            
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
            
            
            %ccSignal{cc,1} = currentBinarySignal;        % binary signal (1 = high, 0 = low)
            
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
            
        end
        clear cc volumes timestamps_hr isDrop curveIDs growthRates_fullOnly 
        clear currentGrowthRates currentVolumes currentTimestamps cellCyles_final 
        clear binaryNutrientSignal nScore currentBinarySignal currentNscore
        clear isSwitch numShifts shiftStage
        
        
        
        % 14. trim cell cycle data to avoid single point cell cycles
        addedVol = ccData(:,3);
        traits = ccData(addedVol > 0,:);
        %signal = ccSignal(addedVol > 0,:);
        clear ccData addedVol ccSignal
        
        
        
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
        switches_final = switches(V_summed == 2,:);
        
        clear birthSize_median birthSize_std divSize_median divSize_std
        clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
        clear V_birth_binary V_division_binary birth_outliers div_outliers V_div V_birth V_summed interdiv
        
        
        
        % plot distribution of volume at division
        
        
        % 16. bin data and normalize bin counts by total counts
        divSize = traits_final(:,2);
        divSize_counts = length(divSize);
        binSize = 0.3;
        assignedBins_raw = ceil(divSize/binSize);
        clear traits
        
        % some conditions do not have data after trimming
        if isempty(traits_final) == 1
            
            filler = linspace(2,6,30);
            
            figure(1)
            plot(filler,zeros(length(filler)),'Color',color,'LineWidth',1)
            hold on
            legend('fluc','low','ave','high')
            title('raw pdf')
            xlabel('division size (cubic um)')
            ylabel('pdf')
            
            
            
            filler = zeros(5,2);
            
            figure(condition+1)
            bar(filler)
            title(strcat('traits of subpopulations: condition-',num2str(condition)))
            
            clear filler color
            
        else
            
            % 17. calculate pdf
            binned_divSize_raw = accumarray(assignedBins_raw, divSize, [], @(x) {x});
            binCounts_raw = cellfun(@length,binned_divSize_raw);
            pdf_divSize_raw = binCounts_raw./divSize_counts;
            clear binned_divSize_raw binCounts_raw divSize_counts
            
            
            % 18. plot PDF
            raw_vals = (1:max(assignedBins_raw)).*binSize; % raw division sizes
            
%             figure(1)
%             plot(raw_vals,pdf_divSize_raw,'Color',color,'LineWidth',1)
%             hold on
%             title('raw pdf')
%             legend('fluc','low','ave','high')
%             xlabel('division size (cubic um)')
%             ylabel('pdf')
            
            clear assignedBins_raw raw_vals pdf_divSize_raw
            
            % split condition data into two populations: above and below average
            
            divSize_mean = mean(divSize);
            divSize_std = std(divSize);
            
            sizeThresh_above = divSize_mean + divSize_std;
            sizeThres_below = divSize_mean - divSize_std;
            
            traits_above = traits_final(divSize > sizeThresh_above,:);
            traits_below = traits_final(divSize < sizeThres_below,:);
            
            
            % compare the mean value for a variety of traits
            
            % for a bar graph with paired entries (one value per subpop):
            % make a matrix with two columns (one for each subpop) and one row per trait
            
            % row 1) interdivision time
            %     2) birth size
            %     3) V_div/V_birth ratio
            %     4) mean growth rate of cell cycle
            %     5) mean growth rate of cell cycle / mean growth rate of population
            %     6) nutrient score
            %     7) cell cycle fraction at the point of a shift (shift stage)
            
           
            
            trait_comparison = [ ...
                
            mean(traits_above(:,6)), mean(traits_below(:,6)); ... % interdiv
            mean(traits_above(:,1)), mean(traits_below(:,1)); ... % birth size
            mean(traits_above(:,2)./traits_above(:,1)), mean(traits_below(:,2)./traits_below(:,1)); ... % V_div / V_birth
            mean(traits_above(:,7)), mean(traits_below(:,7)); ... % mean growth rate
            mean(traits_above(:,9)), mean(traits_below(:,9)); ... % nScore
            nanmean(traits_above(:,10)), nanmean(traits_below(:,10)); ... % numShifts
            nanmean(traits_above(:,11)), nanmean(traits_below(:,11)); ... % cell cycle stage of shift
           
            ];
        
        
            trait_comparison_error = [ ...

            std(traits_above(:,6)), std(traits_below(:,6)); ... % interdiv
            std(traits_above(:,1)), std(traits_below(:,1)); ... % birth size
            std(traits_above(:,2)./traits_above(:,1)), std(traits_below(:,2)./traits_below(:,1)); ... % V_div / V_birth
            std(traits_above(:,7)), std(traits_below(:,7)); ... % mean growth rate
            std(traits_above(:,9)), std(traits_below(:,9)); ... % nScore
            nanstd(traits_above(:,10)), nanstd(traits_below(:,10)); ... % numShifts
            nanstd(traits_above(:,11)), nanstd(traits_below(:,11)); ... % cell cycle stage of shift

            ];
    

            figure(condition+1)
            b = bar(trait_comparison);
            hold on
            errorbar([0.85,1.15; 1.85,2.15; 2.85,3.15; 3.85,4.15; 4.85,5.15; 5.85,6.15; 6.85,7.15],trait_comparison,trait_comparison_error,'.','Color',rgb('Black'))
            title(strcat('traits of subpopulations: condition-',num2str(condition)))
            
            b(1).FaceColor = color;
            b(2).FaceColor = color_b;

            clear trait_comparison trait_comparison_error
    
    
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


