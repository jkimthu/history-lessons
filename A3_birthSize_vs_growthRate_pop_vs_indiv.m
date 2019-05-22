%% Figure A2: mean birth size vs mean growth rate


%  Goal: distinguish single-cell from population-averaged data
%        as in Figure 1C of Taheri et al., Current Biology (2014)

%        single-cell behavior systematically deviates from population-level
%        growth law, relating cell size and growth rate.

%        do flucltuation populations deviate in the same way as
%        steady-state populations?



%  Strategy: 
%
%  Part 1. initialize analysis
%
%       0. initialize experiment data
%       0. define method of calculating growth rate
%       0. initialize colors for plotting
%       0. create array of experiments to use in analysis
%       0. initialize data vectors to store stats for each experiment
%
%  Part 2. collect single cell birth size and instantaneous growth rates
%
%       1.  loop through each experiment to collect data
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. for each condition in experiment...
%                       5. isolate condition specific data

%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. trim condition and growth rate data to include only full cell cycles
%                       9. isolate isDrop, timestamp, and timeSinceBirth data


%       2. inter-division time and mean log2 growth rate of each cycle
%       3. plot population-averaged interdiv time vs mean log2 growth rate



%  Last edit: jen, 2019 May 22
%  Commit: first commit, adapted A2 from history-lessons to include single
%          cell data


%  OK let's go!

%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. create array of experiments to use in analysis
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
compiled_birthSize = cell(length(exptArray),1);
compiled_mu = cell(length(exptArray),1);



%% Part 2. collect single cell birth size and instantaneous growth rates

% 1. loop through each experiment to collect data
for e = 1:length(exptArray)

    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. initialize experiment specific variables
    birthSize = cell(length(bubbletime),1);
    mu_instantaneous = cell(length(bubbletime),1);
    
    
    
    % 4. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    

    
    % for each condition in experiment
    for condition = 1:length(bubbletime)
            
            
        % 5. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end date
        
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes isDrop trackNum timestamps_sec
        
        
        
        % 8. trim condition and growth rate data to include only full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        growthRates_fullOnly = growthRates(curveFinder > 0,:);
        curveIDs_fullOnly = curveFinder(curveFinder > 0,:);   % for trimming growth rate
        curveIDs_unique = unique(curveIDs_fullOnly);          % for assigning birth sizes
        clear conditionData curveFinder growthRates growthRates_all
        
        
        
        % 9. isolate timestamp, isDrop and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        volumes = getGrowthParameter(conditionData_fullOnly,'volume');     % calculated va_vals (cubic um)
        clear timestamps
        
        
        
        % 10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthSize = volumes(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly
        
        
        
        % 11. remove zeros, which occur if no full track data exists at a drop
        Vbirth = final_birthSize(final_birthSize > 0);
        birthTimestamps = finalTimestamps(final_birthSize > 0);
        curveIDs = curveIDs_unique(final_birthSize > 0);
        clear final_birthSize finalTimestamps volumes
        
        
        
        % 12. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            Vbirth_bubbleTrimmed = Vbirth(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            Vbirth_bubbleTrimmed = Vbirth;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear timestamps_hr maxTime isDrop Vbirth curveIDs birthTimestamps
        
        
        
        % 13. truncate data to stabilized regions
        minTime = 3;
        Vbirth_fullyTrimmed = Vbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear Vbirth_bubbleTrimmed curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
        % 14. if no div data in steady-state, skip condition
        if isempty(Vbirth_fullyTrimmed) == 1
            continue
        else
            
            % 14. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            birthSize_median = median(Vbirth_fullyTrimmed);
            birthSize_std_temp = std(Vbirth_fullyTrimmed);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            birthSize_temp = Vbirth_fullyTrimmed(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3)); % cut largest vals, over 3 std out
            IDs_temp = curveIDs_fullyTrimmed_cc(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            birthSize_final = birthSize_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));          % cut smallest vals, over 3 std out
            IDs_final = IDs_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));   
            clear birthSize_median birthSize_std_temp birthSize_temp IDs_temp
            
            % iv. remove corresponding growth rates from datasets
            trimmedIDs = setdiff(curveIDs_unique,IDs_final);    % curve IDs in growth rate dataset, NOT in final IDs trimmed by cell cycle
            toTrim = ismember(curveIDs_fullOnly,trimmedIDs);   % vector of what to trim or not in growth rate
            trimmed_curves_insta = curveIDs_fullOnly(toTrim == 0);
            trimmed_mus = growthRates_fullOnly(toTrim == 0);
            clear toTrim trimmedIDs curveIDs_fullOnly growthRates_fullOnly
            
            
            
            % 15. bin growth rates by cell cycle, to match organization of birth size data
            mus_binned = accumarray(trimmed_curves_insta,trimmed_mus,[],@(x) {x});
            mus = mus_binned(~cellfun('isempty',mus_binned));
            clear trimmed_curves_insta trimmed_mus
            
            
            
            % 16. store condition data into one variable per experiment
            birthSize{condition} = birthSize_final;
            mu_instantaneous{condition} = mus;
        
            
        end
        
    end
      
    % 17. store experiment data into single variable for further analysis
    compiled_birthSize{e} = birthSize;
    compiled_mu{e} = mu_instantaneous;
    
end

% 18. save hard earned data
save('A3_data.mat','compiled_birthSize','compiled_mu','exptArray')



%% Part 3. plotting

% 0. initialize plotting parameters
xmin = -1;                      % lower limit for plotting x axis
xmax = 5;                       % upper limit for plotting x axis

palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
shapes = {'o','x','square','*'};



%%

rando lines for plotting PDFs

 
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
        errorbar(gr_mean{e},Vbirth_mean{e},Vbirth_std{e},'Color',color)
        hold on
        plot(gr_mean{e},Vbirth_mean{e},'Marker',xmark,'Color',color)
        hold on
        ylabel('birth volume (cubic um)')
        xlabel('growth rate (1/hr)')
        title('population mean from all experiments')
        axis([0 4 0 8])

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