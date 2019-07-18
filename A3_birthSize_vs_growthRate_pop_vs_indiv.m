%% Figure A3: mean birth size vs mean growth rate


%  Goal: distinguish single-cell from population-averaged data
%        as in Figure 1C of Taheri et al., Current Biology (2014)

%        single-cell behavior systematically deviates from population-level
%        growth law, relating cell size and growth rate.

%        do flucltuation populations deviate in the same way as
%        steady-state populations?



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell birth size and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line




%  Last edit: jen, 2019 July 18
%  Commit: edit to plot best fit line over steady points and calculate CV


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
compiled_birthLength = cell(length(exptArray),1);
compiled_mu = cell(length(exptArray),1);



%% Part 2. collect single cell birth size and instantaneous growth rates


%  Strategy:
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
%                      10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
%                      11. remove zeros, which occur if no full track data exists at a drop
%                      12. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      13. truncate data to stabilized regions
%                      14. if no div data in steady-state, skip condition
%                          else, trim outliers (those 3 std dev away from median) from final dataset
%                      15. bin growth rates by cell cycle, to match organization of birth size data
%                      16. store condition data into one variable per experiment
%               17. store experiment data into single variable for further analysis
%      18. save hard earned data


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
        
        
        
        % 9. isolate timestamp, isDrop, length and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        volumes = getGrowthParameter(conditionData_fullOnly,'volume');     % calculated va_vals (cubic um)
        majorAxis = getGrowthParameter(conditionData_fullOnly,'length');   % length (um)
        clear timestamps
        
        
        
        % 10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthSize = volumes(isDrop==1);
        final_birthLength = majorAxis(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly
        
        
        
        % 11. remove zeros, which occur if no full track data exists at a drop
        Vbirth = final_birthSize(final_birthSize > 0);
        Lbirth = final_birthLength(final_birthSize > 0);
        birthTimestamps = finalTimestamps(final_birthSize > 0);
        curveIDs = curveIDs_unique(final_birthSize > 0);
        clear final_birthSize finalTimestamps volumes majorAxis
        
        
        
        % 12. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            Vbirth_bubbleTrimmed = Vbirth(birthTimestamps <= maxTime,:);
            Lbirth_bubbleTrimmed = Lbirth(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            Vbirth_bubbleTrimmed = Vbirth;
            Lbirth_bubbleTrimmed = Lbirth;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear timestamps_hr maxTime isDrop Vbirth Lbirth curveIDs birthTimestamps
        
        
        
        % 13. truncate data to stabilized regions
        minTime = 3;
        Vbirth_fullyTrimmed = Vbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        Lbirth_fullyTrimmed = Lbirth_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear Vbirth_bubbleTrimmed Lbirth_bubbleTrimmed curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
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
            birthLength_temp = Lbirth_fullyTrimmed(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3));
            IDs_temp = curveIDs_fullyTrimmed_cc(Vbirth_fullyTrimmed <= (birthSize_median+birthSize_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            birthSize_final = birthSize_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));          % cut smallest vals, over 3 std out
            birthLength_final = birthLength_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3)); 
            IDs_final = IDs_temp(birthSize_temp >= (birthSize_median-birthSize_std_temp*3));   
            clear birthSize_median birthSize_std_temp birthSize_temp birthLength_temp IDs_temp
            
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
            birthLength{condition} = birthLength_final;
            mu_instantaneous{condition} = mus;
        
            
        end
        
    end
      
    % 17. store experiment data into single variable for further analysis
    compiled_birthSize{e} = birthSize;
    compiled_birthLength{e} = birthLength;
    compiled_mu{e} = mu_instantaneous;
    
end

% 18. save hard earned data
save('A3_data_2.mat','compiled_birthSize','compiled_birthLength','compiled_mu','exptArray')



%% Part 3. plotting

% goal: plot of newborn cell size vs growth rate,
%       overlaying population and single cell data

% strategy: 
%
%       i. accumulate data from each condition
%          conditions: each steady and each fluctuating timescale (7 total)
%
%      ii. calculate mean birth size and growth rate for each condition (population data)
%          plot each condition as a closed orange point and fit a line (the growth law ACROSS conditions)
%
%     iii. bin individual data by growth rate
%
%      iv. calculate mean birth size and growth rate for each bin (individual data)
%          plot each bin as an open blue point and fit a line WITHIN each condition


clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A3_data_2.mat')

% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};

shape = 'o';


% 1. accumulate data from each condition
fluc = 1; % row number in data structure
low = 2; 
ave = 3; 
high = 4;

sigmas = 1;

for ee = 1:length(environment_order)
    
    condition = environment_order{ee};
    
    if ischar(condition) == 1
        
        % steady environment! concatenate data based on nutrient level
        if strcmp(condition,'low') == 1
            
            mu_low = [];
            birthSize_low = [];
            birthLength_low = [];
            
            % loop through all experiments and store low data
            for expt = 1:length(compiled_mu)
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{low,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{low,1};
                expt_length = compiled_birthLength{expt,1}{1,low};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{low,1});
                end
                
                % concanetate individual cell cycle values
                mu_low = [mu_low; expt_mu];
                birthSize_low = [birthSize_low; expt_size];
                birthLength_low = [birthLength_low; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_low;
                condition_size = birthSize_low;
                condition_length = birthLength_low;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std * sigmas);
                upper = condition_mu > condition_mean - (condition_std * sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                size{1} = range_size;
                mu{1} = range_mu;
                bLength{1} = range_length;
                
        elseif strcmp(condition,'ave') == 1
            
            mu_ave = [];
            birthSize_ave = [];
            birthLength_ave = [];
            
            % loop through all experiments and store ave data
            for expt = 1:length(compiled_mu)
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{ave,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{ave,1};
                expt_length = compiled_birthLength{expt,1}{1,ave};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{ave,1});
                end
                
                % concanetate individual cell cycle values
                mu_ave = [mu_ave; expt_mu];
                birthSize_ave = [birthSize_ave; expt_size];
                birthLength_ave = [birthLength_ave; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_ave;
                condition_size = birthSize_ave;
                condition_length = birthLength_ave;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std * sigmas);
                upper = condition_mu > condition_mean - (condition_std * sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{6} = range_mu;
                size{6} = range_size;
                bLength{6} = range_length;
            
        elseif strcmp(condition,'high') == 1
            
            mu_high = [];
            birthSize_high = [];
            birthLength_high = [];
            
            % loop through all experiments and store high data
            for expt = 1:length(compiled_mu)
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{high,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{high,1};
                expt_length = compiled_birthLength{expt,1}{1,high};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{high,1});
                end
                
                % concanetate individual cell cycle values
                mu_high = [mu_high; expt_mu];
                birthSize_high = [birthSize_high; expt_size];
                birthLength_high = [birthLength_high; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_high;
                condition_size = birthSize_high;
                condition_length = birthLength_high;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std * sigmas);
                upper = condition_mu > condition_mean - (condition_std * sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{7} = range_mu;
                size{7} = range_size;
                bLength{7} = range_length;
            
        end
    else
        
        % fluctuating environment! concatenate based on timescale
        if condition == 30
            idx = [2,3,4]; % ID of experiments with this fluc timescale
            
            mu_30 = [];
            birthSize_30 = [];
            length_30 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{fluc,1};
                expt_length = compiled_birthLength{expt,1}{1,fluc};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
                end
                
                % concanetate individual cell cycle values
                mu_30 = [mu_30; expt_mu];
                birthSize_30 = [birthSize_30; expt_size];
                length_30 = [length_30; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_30;
                condition_size = birthSize_30;
                condition_length = length_30;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std * sigmas);
                upper = condition_mu > condition_mean - (condition_std * sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{2} = range_mu;
                size{2} = range_size;
                bLength{2} = range_length;
            
            
        elseif condition == 300 
            idx = [5,6,7]; % ID of experimennts with this fluc timescale
            
            mu_300 = [];
            birthSize_300 = [];
            length_300 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{fluc,1};
                expt_length = compiled_birthLength{expt,1}{1,fluc};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
                end
                
                % concanetate individual cell cycle values
                mu_300 = [mu_300; expt_mu];
                birthSize_300 = [birthSize_300; expt_size];
                length_300 = [length_300; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_300;
                condition_size = birthSize_300;
                condition_length = length_300;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std * sigmas);
                upper = condition_mu > condition_mean - (condition_std * sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{3} = range_mu;
                size{3} = range_size;
                bLength{3} = range_length;
            
                
        elseif condition == 900
            idx = [9,10,11,12]; % ID of experimennts with this fluc timescale
            
            mu_900 = [];
            birthSize_900 = [];
            length_900 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{fluc,1};
                expt_length = compiled_birthLength{expt,1}{1,fluc};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
                end
                
                % concanetate individual cell cycle values
                mu_900 = [mu_900; expt_mu];
                birthSize_900 = [birthSize_900; expt_size];
                length_900 = [length_900; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_900;
                condition_size = birthSize_900;
                condition_length = length_900;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std*sigmas);
                upper = condition_mu > condition_mean - (condition_std*sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{4} = range_mu;
                size{4} = range_size;
                bLength{4} = range_length;
                
            
        elseif condition == 3600
            idx = [13,14,15]; % ID of experimennts with this fluc timescale
            
            mu_3600 = [];
            birthSize_3600 = [];
            length_3600 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                
                % isolate data 
                expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
                expt_size = compiled_birthSize{expt,1}{fluc,1};
                expt_length = compiled_birthLength{expt,1}{1,fluc};
                
                % if data exists, calculate mean of each cell cycle
                if ~isempty(expt_mu)
                    expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
                end
                
                % concanetate individual cell cycle values
                mu_3600 = [mu_3600; expt_mu];
                birthSize_3600 = [birthSize_3600; expt_size];
                length_3600 = [length_3600; expt_length];
                clear expt_mu expt_size expt_length
                
            end
            
                condition_mu = mu_3600;
                condition_size = birthSize_3600;
                condition_length = length_3600;
                
                % isolate cycles within 3 st dev of mean
                condition_mean = nanmean(condition_mu);
                condition_std = nanstd(condition_mu);
                
                lower = condition_mu < condition_mean + (condition_std*sigmas);
                upper = condition_mu > condition_mean - (condition_std*sigmas);
                combined = lower + upper;
                
                range_mu = condition_mu(combined == 2);
                range_size = condition_size(combined == 2);
                range_length = condition_length(combined == 2);
                clear condition_mean condition_std lower upper
                
                % store condition data
                mu{5} = range_mu;
                size{5} = range_size;
                bLength{5} = range_length;
            
        end  
    end
end
clear fluc low ave high idx condition ee expt arrayIndex
clear condition_mu condition_size condition_length
clear range_mu range_size range_length




% 2. calculate mean birth size and growth rate for each condition (population data)
%    plot each condition as a closed orange point and fit a line (the growth law ACROSS conditions)


% assemble data based on environment order
population_birthSize = cellfun(@mean,size); % population birth size of each condition
population_mu = cellfun(@nanmean,mu);       % population growth rate of each condition
population_birthLength = cellfun(@mean,bLength); % population birth length of each condition




% plot each population point
figure(1)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthSize(cc)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
    hold on
    
end
clear color cc

figure(2)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthLength(cc)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
    hold on
    
end
clear color cc




% 3. bin individual data by growth rate
binFactor = 10; % each bin is 0.25 1/h

for condition = 1:length(environment_order)
    
    condition_color = rgb(palette(condition));
    condition_mu = mu{condition};
    condition_size = size{condition};
    condition_length = bLength{condition};
    
    bins = floor((condition_mu + 10)* binFactor);
    binned_mu = accumarray(bins,condition_mu,[],@mean);
    indiv_mu = binned_mu(binned_mu > 0);
    
    binned_size = accumarray(bins,condition_size,[],@mean);
    indiv_size = binned_size(binned_mu > 0);
    binned_length = accumarray(bins,condition_length,[],@mean);
    indiv_length = binned_length(binned_mu > 0);
    
    
    % 4. calculate mean birth size and growth rate for each bin (individual data)
    % plot each bin as an open blue point and fit a line WITHIN each condition
    figure(1)
    hold on
    plot(indiv_mu,log(indiv_size),'Color',condition_color,'Marker',shape,'MarkerSize',6,'LineWidth',1)
    
    
    figure(2)
    hold on
    plot(indiv_mu,log(indiv_length),'Color',condition_color,'Marker',shape,'MarkerSize',6,'LineWidth',1)
    
end
clear condition_color condition_mu condition_size condition_length
clear bins binned_mu indiv_mu binned_size indiv_size binned_length indiv_length


figure(1)
legend('low','30','300','900','3600','ave','high')
title('population vs individual growth laws')
xlabel('mean mu')
ylabel('mean birth volume')
ylim([0.5 2])
xlim([-0.1 7])

figure(2)
legend('low','30','300','900','3600','ave','high')
title('population vs individual growth laws')
xlabel('mean mu')
ylabel('mean birth length')
ylim([0.5 1.75])
xlim([-0.1 7])


% save hard earned data!
save('A3_data_1sigma.mat','compiled_birthSize','compiled_birthLength','compiled_mu','exptArray','population_birthSize','population_mu','population_birthLength','size','mu','bLength')



%% Part 4. fitting 

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A3_data_2.mat') % as of 2019-07-18, saved data compiles cc within 3 sigmas

palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};




% 1. determine which data points are steady = 1 or fluctuating = 0
for eo = 1:length(environment_order)
    isSteady(eo) = ischar(environment_order{eo});
end
clear eo


% 2. isolate steady data
steady_mu = population_mu(isSteady == 1);
steady_birthLength = population_birthLength(isSteady == 1);
steady_birthVolume = population_birthSize(isSteady == 1);


% 3. best fit line
fit_length = polyfit(steady_mu,log(steady_birthLength),1);
fit_vol = polyfit(steady_mu,log(steady_birthVolume),1);

% 4. plot data
figure(1)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthSize(cc)),'Color',color,'Marker','o','MarkerSize',10,'LineWidth',2)
    hold on
    
end

figure(2)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthLength(cc)),'Color',color,'Marker','o','MarkerSize',10,'LineWidth',2)
    hold on
    
end
clear color cc


% 5. plot best fit line
x = linspace(steady_mu(1),steady_mu(end),10);
yv = fit_vol(1)*x + fit_vol(2);
yl = fit_length(1)*x + fit_length(2);

figure(1)
hold on
plot(x,yv,'Color',rgb('SlateGray'))
axis([0.5 3.5 0.6 1.8])
xlabel('mean mu')
ylabel('ln(mean birth volume)')
legend('low','30','300','900','3600','ave','high')
text(2, 1.65, strcat('y=',num2str(fit_vol(1)),'x+',num2str(fit_vol(2))))

figure(2)
hold on
plot(x,yl,'Color',rgb('SlateGray'))
axis([0.5 3.5 0.6 1.6])
xlabel('mean mu')
ylabel('ln(mean birth length)')
legend('low','30','300','900','3600','ave','high')
text(2, 1.35, strcat('y=',num2str(fit_length(1)),'x+',num2str(fit_length(2))))



%% 5. coefficient of variation analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')


% 0. determine number of sigmas in dataset
%sigs = '3sigmas';
sigs = '1sigma';


% 0. load dataset
load(strcat('A3_data_',sigs,'.mat')) % as of 2019-07-18, saved data compiles cc within 3 sigmas



% 1. for 3 sigma dataset, calculate cv for each condition
pop_std_size = cellfun(@std,size);
pop_std_mu = cellfun(@nanstd,mu);
pop_std_birthLength = cellfun(@std, bLength);


% 2. calculate CV for each condition
cv_volume_3sigma = pop_std_size./population_birthSize;
cv_length_3sigma = pop_std_birthLength./population_birthLength;
cv_lambda_3sigma = pop_std_mu./population_mu;


% 3. plot CVs
figure(1)
bar([cv_volume_3sigma;cv_length_3sigma;cv_lambda_3sigma])
ylabel('coefficient of variation')
name = {'birth volume';'birth length';'growth rate'};
set(gca,'xticklabel',name)
legend({'low','30 s','5 min','15 min','60 min','ave','high'})
title(sigs)