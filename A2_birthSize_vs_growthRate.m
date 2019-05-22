%% Figure A2: mean birth size vs mean growth rate


%  Goal: determine birth size as predicted
%        by population-averaged growth rate d(logV)/dt, plotting:


%  Strategy: 
%       0. initialize experiment data
%       1. create array of experiments to use in analysis, then loop through each
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. initialize colors for plotting, then loop through conditions
%                       6. isolate condition specific data
%                       7. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       8. calculate growth rate
%                       9. trim condition and growth rate data to include only full cell cycles
%                      10. isolate isDrop, timestamp, and timeSinceBirth data
%                      11. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
%                      12. remove zeros, which occur if no full track data exists at a drop
%                      13. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      14. truncate data to stabilized regions
%                      15. trim outliers (those 3 std dev away from median) from final dataset
%                      16. calculate final populuation mean and count
%                      17. plot mean birth size vs mean growth rate




%  Last edit: jen, 2019 May 22
%  Commit: minor tidying while working on A3


%  OK let's go!

%% initialize

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates
xmin = -1;                      % lower limit for plotting x axis
xmax = 5;                       % upper limit for plotting x axis


%%

% 1a. create array of experiments to use in analysis, then loop through each
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 1b. initialize data vectors to store stats for each experiment
Vbirth_mean = cell(length(exptArray),1);
Vbirth_std = cell(length(exptArray),1);
Vbirth_count = cell(length(exptArray),1);

gr_mean = cell(length(exptArray),1);
gr_std = cell(length(exptArray),1);
gr_count = cell(length(exptArray),1);


for e = 1:length(exptArray)

    
    % 2. initialize experiment meta data
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
    clear D5 T xy_start xy_end expType date
    
    
    % 5. initialize colors for plotting, then loop through conditions
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    
    
    for condition = 1:length(bubbletime)
        
        
        % 6. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        
        
        
        % 7. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
        
        
        
        
        % 8. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes timestamps_sec isDrop curveFinder trackNum
        
        
        
        % 9. trim condition and growth rate data to include only full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        growthRates_fullOnly = growthRates(curveIDs > 0,:);
        clear curveFinder conditionData curveIDs growthRates
        
        
        
        % 10. isolate isDrop, timestamp, and timeSinceBirth data
        timestamps_hr = conditionData_fullOnly(:,2)./3600;    % col 2  = timestamp in seconds converted to hours
        timestamps_hr2 = conditionData_fullOnly(:,2)./3600;   % same as above, but for growth rate trimming
        
        isDrop = conditionData_fullOnly(:,4);                 % col 4  = isDrop, 1 marks a birth event
        volumes = conditionData_fullOnly(:,11);             % col 11 = calculated va_vals (cubic um)
        
        
        
        % 11. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthSize = volumes(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        
        
        
        % 12. remove zeros, which occur if no full track data exists at a drop
        Vbirth = final_birthSize(final_birthSize>0);
        birthTimestamps = finalTimestamps(final_birthSize>0);
        clear final_birthSize finalTimestamps timeSinceBirth_min
        
        
        
        % 13. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        
        if maxTime > 0
            Vbirth_bubbleTrimmed = Vbirth(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
            
            growthRates_bubbleTrimmed = growthRates_fullOnly(timestamps_hr2 <= maxTime,:);
            grTimestamps_bubbleTrimmed = timestamps_hr2(timestamps_hr2 <= maxTime,:);
            
        else
            Vbirth_bubbleTrimmed = Vbirth;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
            growthRates_bubbleTrimmed = growthRates_fullOnly;
            grTimestamps_bubbleTrimmed = timestamps_hr2;
        end
        clear timestamps_hr timestamps_hr2 maxTime isDrop birthIndeces divTimestamps divTimes growthRates_fullOnly
        
        
        
        % 14. truncate data to stabilized regions
        minTime = 3;
        birthSize_fullyTrimmed = Vbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        growthRates_fullyTrimmed = growthRates_bubbleTrimmed(grTimestamps_bubbleTrimmed >= minTime,:);
        clear growthRates_bubbleTrimmed divTimes_bubbleTrimmed divTimestamps_bubbleTrimmed grTimestamps_bubbleTrimmed
        
        
        
        % 15. if no div data in steady-state, skip condition
        if isempty(birthSize_fullyTrimmed) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            birthSize_median = median(birthSize_fullyTrimmed);
            birthSize_std_temp = std(birthSize_fullyTrimmed);
            birthSize_temp = birthSize_fullyTrimmed(birthSize_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3)); % cut smallest vals, over 3 std out
            birthSize_final = birthSize_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));          % cut largest vals, over 3 std out
            clear birthSize_median birthSize_std_temp birthSize_temp
            
            
            % 16. calculate final populuation mean and count
            Vbirth_mean{e} = mean(birthSize_final);
            Vbirth_std{e} = std(birthSize_final);
            Vbirth_count{e} = length(birthSize_final);
            
            gr_mean{e} = nanmean(growthRates_fullyTrimmed);
            gr_std{e} = nanstd(growthRates_fullyTrimmed);
            gr_count{e} = sum(isnan(growthRates_fullyTrimmed)==0);
            
        end
        
        
        % 17. plot mean birth size vs mean growth rate
        color = rgb(palette(condition));
        
        if condition == 1 && timescale == 300
            xmark = shapes{2};
        elseif condition == 1 && timescale == 900
            xmark = shapes{3};
        elseif condition == 1 && timescale == 3600
            xmark = shapes{4};
        else
            xmark = shapes{1};
        end
        
        
        figure(1)
        errorbar(gr_mean{e},Vbirth_mean{e},Vbirth_std{e},'Color',color)
        hold on
        plot(gr_mean{e},Vbirth_mean{e},'Marker',xmark,'Color',color)
        hold on
        ylabel('birth volume (cubic um)')
        xlabel('growth rate (1/hr)')
        title('population mean from all experiments')
        axis([0 4 0 8])

        
    end
      
    
end


%% rando lines for plotting PDFs

% 16. bin data and normalize bin counts by total counts
assignedBins_raw = ceil(birthSize_final/binSize);



% 15. calculate pdf
binned_divT_raw = accumarray(assignedBins_raw, birthSize_final, [], @(x) {x});
binCounts_raw = cellfun(@length,binned_divT_raw);
pdf_divT_raw = binCounts_raw/Vbirth_count;



% 16. plot PDF
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
color = rgb(palette(condition));

% raw growth rates
raw_vals = (1:max(assignedBins_raw)).*binSize;

figure(1)
plot(raw_vals,pdf_divT_raw,'Color',color,'LineWidth',1)
hold on
title('raw pdf')
legend('fluc','low','ave','high')
xlabel('interdivision time (min)')
ylabel('pdf')