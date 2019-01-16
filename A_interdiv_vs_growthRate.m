%% Figure A: mean interdivision time vs mean growth rate


%  Goal: determine interdivision time as predicted
%        by population-averaged growth rate d(logV)/dt, plotting:


%  Strategy: 
%       0. initialize experiment data
%       1. create array of experiments to use in analysis, then loop through each
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. initialize colors for plotting, then loop through conditions
%                       5. isolate condition specific data
%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. trim condition and growth rate data to include only full cell cycles
%                       9. isolate isDrop, timestamp, and timeSinceBirth data


%       2. inter-division time and mean log2 growth rate of each cycle
%       3. plot population-averaged interdiv time vs mean log2 growth rate



%  Last edit: jen, 2019 January 16
%  Commit: first commit, adapted from figure 9 from sixpack-abs


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
divT_mean = cell(length(exptArray),1);
divT_std = cell(length(exptArray),1);
divT_count = cell(length(exptArray),1);

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
    if strcmp(date,'2017-11-12') == 1
        filename = strcat('lb-fluc-',date,'-width1p4-jiggle-0p5.mat');
    elseif strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','T');
    
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType date
    
    
    % 4. initialize colors for plotting, then loop through conditions
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    
    
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = conditionData(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes timestamps_sec isDrop curveFinder trackNum
        
        
        
        % 8. trim condition and growth rate data to include only full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        growthRates_fullOnly = growthRates(curveIDs > 0,:);
        clear curveFinder conditionData curveIDs growthRates
        
        
        
        % 9. isolate isDrop, timestamp, and timeSinceBirth data
        timestamps_hr = conditionData_fullOnly(:,2)./3600;    % col 2  = timestamp in seconds converted to hours
        timestamps_hr2 = conditionData_fullOnly(:,2)./3600;   % same as above, but for growth rate trimming
        
        isDrop = conditionData_fullOnly(:,4);                 % col 4  = isDrop, 1 marks a birth event
        timeSinceBirth_min = conditionData_fullOnly(:,6)./60; % col 6  = timeSinceBirth in seconds converted to min
        
        
        
        % 10. extract only final timeSinceBirth from each growth curve,
        %     this is the inter-division time!
        birthIndeces = find(isDrop==1);
        finalDivTimes = timeSinceBirth_min(birthIndeces(2:end)-1); % final timeSince birth is indexed right before every drop event
               % ^ minutes between recorded birth and division
               
        finalTimestamps = timestamps_hr(birthIndeces(2:end)-1); % experiment timestamp (hours) of each division event.
        
        
        
        % 11. remove zeros, which occur if no full track data exists at a drop
        divTimes = finalDivTimes(finalDivTimes>0);
        divTimestamps = finalTimestamps(finalDivTimes>0);
        clear finalDivTimes finalTimestamps birthIndeces timeSinceBirth_min
        
        
        
        % 12. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        
        if maxTime > 0
            divTimes_bubbleTrimmed = divTimes(divTimestamps <= maxTime,:);
            divTimestamps_bubbleTrimmed = divTimestamps(divTimestamps <= maxTime,:);
            
            growthRates_bubbleTrimmed = growthRates_fullOnly(timestamps_hr2 <= maxTime,:);
            grTimestamps_bubbleTrimmed = timestamps_hr2(timestamps_hr2 <= maxTime,:);
            
        else
            divTimes_bubbleTrimmed = divTimes;
            divTimestamps_bubbleTrimmed = divTimestamps;
            
            growthRates_bubbleTrimmed = growthRates_fullOnly;
            grTimestamps_bubbleTrimmed = timestamps_hr2;
        end
        clear timestamps_hr timestamps_hr2 maxTime isDrop birthIndeces divTimestamps divTimes growthRates_fullOnly
        
        
        
        % 13. truncate data to stabilized regions
        minTime = 3;
        divTimes_fullyTrimmed = divTimes_bubbleTrimmed(divTimestamps_bubbleTrimmed >= minTime,:);
        growthRates_fullyTrimmed = growthRates_bubbleTrimmed(grTimestamps_bubbleTrimmed >= minTime,:);
        clear growthRates_bubbleTrimmed divTimes_bubbleTrimmed divTimestamps_bubbleTrimmed grTimestamps_bubbleTrimmed
        
        
        
        % if no div data in steady-state, skip condition
        if isempty(divTimes_fullyTrimmed) == 1
            continue
        else
            
            % 14. trim outliers (those 3 std dev away from median) from final dataset
            divT_median = median(divTimes_fullyTrimmed);
            divT_std_temp = std(divTimes_fullyTrimmed);
            divT_temp = divTimes_fullyTrimmed(divTimes_fullyTrimmed <= (divT_median+divT_std_temp*3)); % cut smallest vals, over 3 std out
            divT_final = divT_temp(divT_temp >= (divT_median-divT_std_temp*3));          % cut largest vals, over 3 std out
            clear divT_median divT_std_temp divT_temp
            
            
            % 15. calculate final populuation mean and count
            divT_mean{e} = mean(divT_final);
            divT_std{e} = std(divT_final);
            divT_count{e} = length(divT_final);
            
            gr_mean{e} = nanmean(growthRates_fullyTrimmed);
            gr_std{e} = nanstd(growthRates_fullyTrimmed);
            gr_count{e} = sum(isnan(growthRates_fullyTrimmed)==0);
            
        end
        
        
        % 16. plot mean interdivT vs mean growth rate
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
        errorbar(gr_mean{e},divT_mean{e},divT_std{e},'Color',color)
        hold on
        plot(gr_mean{e},divT_mean{e},'Marker',xmark,'Color',color)
        hold on
        ylabel('inter-division time (min)')
        xlabel('growth rate (1/hr)')
        title('population mean from all experiments')
        axis([0 19 0 80])

        
    end
      
    
end


%% rando lines for plotting PDFs

% 16. bin data and normalize bin counts by total counts
assignedBins_raw = ceil(divT_final/binSize);



% 15. calculate pdf
binned_divT_raw = accumarray(assignedBins_raw, divT_final, [], @(x) {x});
binCounts_raw = cellfun(@length,binned_divT_raw);
pdf_divT_raw = binCounts_raw/divT_count;



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