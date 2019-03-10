%% Figure B1: mean V_div  vs mean V_birth


%  Goal: plot empirical relationship bewteen size at division and at birth.
%        plot population-averaged data from each replicate condition.


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) trim data to post-stabilized and pre-bubble time
%       d) trim data to only consider curves within 3 std of mean interdiv times
%       e) collect V_division and V_birth



%  Last edit: jen, 2019 Feb 27

%  Commit: first commit, checking 2x relationship between Ldiv and Lbirth


%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
Ldiv_mean = cell(length(exptArray),1);
Ldiv_std = cell(length(exptArray),1);
Vdiv_count = cell(length(exptArray),1);

Lbirth_mean = cell(length(exptArray),1);
Lbirth_std = cell(length(exptArray),1);
Vbirth_count = cell(length(exptArray),1);


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
shapes = {'o','x','square','*'};

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
    
    datesForLegend{e} = date;
    
    
    % 3. load measured experiment data    
    %experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    %cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    % 4. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType date
    
    
    
    % for each condition in experiment...
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        
                
        % 6. trim data to full cell cycles ONLY
        curveFinder = getGrowthParameter(conditionData,'curveFinder');     % curveFinder = ID num of cell cycles bound by divisions
        fullData = conditionData(curveFinder > 0,:);
        clear curveFinder data_fullyTrimmed data_bubbleTrimmed
        
        
        % 7. isolate volume, isDrop, curveID and timestamp data
        lengths = getGrowthParameter(fullData,'length');
        isDrop = getGrowthParameter(fullData,'isDrop');     % isDrop == 1 marks a birth event
        curveID = getGrowthParameter(fullData,'curveFinder');
        timestamps_sec = getGrowthParameter(fullData,'timestamp');
        timestamps_hr = timestamps_sec./3600;               % timestamp in seconds converted to hours, for time-based trim
        clear timestamps_sec
        
        
        % 8. identity unique cell cycles by ID number
        cellCycles = curveID(isDrop == 1);
        birthTimes = timestamps_hr(isDrop == 1);
        
        
        % 9. remove birth times prior to 3 hr
        birthTimes_post3 = birthTimes(birthTimes > 3); 
        cellCycles_post3 = cellCycles(birthTimes > 3); 
        clear cellCycles birthTimes
        
        
        % 10. remove birth times post bubbles
        birthTimes_final = birthTimes_post3(birthTimes_post3 < bubbletime(condition));
        cellCyles_final = cellCycles_post3(birthTimes_post3 < bubbletime(condition));
        clear cellCycles_post3 birthTimes_post3 birthTimes_final
        
        
        % 11. for remaining cell cycles, identify volume at birth and at division of each
        ccData = nan(length(cellCyles_final),6);
        
        for cc = 1:length(cellCyles_final)
            
            % isolate volume and timestamps for current cell cycle
            currentLengths = lengths(curveID == cellCyles_final(cc));
            currentTimestamps = timestamps_hr(curveID == cellCyles_final(cc));
            
            ccData(cc,1) = currentLengths(1);       % V_birth
            ccData(cc,2) = currentLengths(end);     % V_division
            ccData(cc,3) = currentLengths(end) - currentLengths(1); % added volume
            ccData(cc,4) = currentTimestamps(1);    % T_birth
            ccData(cc,5) = currentTimestamps(end);  % T_division
            ccData(cc,6) = currentTimestamps(end) - currentTimestamps(1); % interdivision time
            
        end
        clear cc lengths timestamps_hr isDrop curveID currentLengths currentTimestamps cellCyles_final
        
        
        % 12. isolate inter-division time data and calculate stats
        interdivTime_hr = ccData(:,6);
        interdivTime_hr(interdivTime_hr == 0) = NaN;       % some zeros creep in, remove
        divTime_median = nanmedian(interdivTime_hr);
        divTime_std = nanstd(interdivTime_hr);

        
        
        % 12. trim cell cycles outside of 3 st dev of median
        upper = divTime_median + (3*divTime_std);
        lower = divTime_median - (3*divTime_std);
        ccData_temp = ccData((interdivTime_hr < lower == 0),:);
        ccData_final = ccData_temp((ccData_temp(:,6) > upper == 0),:);
        clear ccData_temp upper lower divTime_median divTime_std interdivTime_hr
        

        
        % 13. trim outliers (those 3 std dev away from median) from final dataset
        %     a) id median and std of division size
        %     b) id median and std of birth size
        %     c) find indeces in both vectors that are within 3 std
        %     d) keep data from indeces in both vectors
        L_birth = ccData_final(:,1);
        L_div = ccData_final(:,2);
        interdiv = ccData_final(:,6);
        
        divSize_median = median(L_div); 
        divSize_std = std(L_div);
        
        birthSize_median = median(L_birth); 
        birthSize_std = std(L_birth);
        
        div_bigOutlier = find(L_div > (divSize_median+divSize_std*3));
        div_smallOutlier = find(L_div < (divSize_median-divSize_std*3));
        div_outliers = [div_bigOutlier; div_smallOutlier];
        
        birth_bigOutlier = find(L_birth > (birthSize_median+birthSize_std*3));
        birth_smallOutlier = find(L_birth < (birthSize_median-birthSize_std*3));
        birth_outliers = [birth_bigOutlier; birth_smallOutlier];

        L_division_binary = ones(length(L_div),1);
        L_division_binary(div_outliers) = 0;
        
        L_birth_binary = ones(length(L_birth),1);
        L_birth_binary(birth_outliers) = 0;
        
        L_summed = L_division_binary + L_birth_binary;
        
        L_division_final = L_div(L_summed == 2);
        L_birth_final = L_birth(L_summed == 2);
        interdiv_final = interdiv(L_summed == 2);
       
        clear birthSize_median birthSize_std divSize_median divSize_std
        clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
        clear L_birth_binary L_division_binary birth_outliers div_outliers L_div L_birth L_summed interdiv
        
        
                
        % 13. calculate mean and standard deviation of birth and division volume
        Lbirth_mean{e,condition} = mean(L_birth_final);
        Lbirth_std{e,condition} = std(L_birth_final);
        counts{3} = length(L_birth_final);
        
        Ldiv_mean{e,condition} = mean(L_division_final);
        Ldiv_std{e,condition} = std(L_division_final);
        
        interdiv_final_min = interdiv_final.*60;
        interdiv_mean{e,condition} = mean(interdiv_final_min);
        interdiv_std{e,condition} = std(interdiv_final_min);
        clear interdiv_final
        
        
        
        % 14. initialize plotting parameters
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
        
      
     
        
        % 15. plot!
        % (i) division size vs. birth size
        figure(1)
        errorbar(Lbirth_mean{e,condition},Ldiv_mean{e,condition},Ldiv_std{e,condition},'Color',color)
        hold on
        plot(Lbirth_mean{e,condition},Ldiv_mean{e,condition},'Marker',xmark,'Color',color)
        hold on
        xlabel('birth length (um)')
        ylabel('division length (cubic um)')
        title('population averages from all experiments')
        axis([1 5 0 12])
        
        
        % (ii) inter-division time vs. birth size
        figure(2)
        errorbar(Lbirth_mean{e,condition},interdiv_mean{e,condition},interdiv_std{e,condition},'Color',color)
        hold on
        plot(Lbirth_mean{e,condition},interdiv_mean{e,condition},'Marker',xmark,'Color',color)
        hold on
        xlabel('birth length (um)')
        ylabel('inter-division time (min)')
        title('population averages from all experiments')
        axis([1 5 0 80])
        
        
    end
    
    
end


% expectations
expected_color = rgb('Silver');
expected_Size_b = linspace(1,8,8);
expected_Sive_div = 2* expected_Size_b ;

figure(1)
hold on
plot(expected_Size_b,expected_Sive_div,'-','Color',expected_color);
saveas(gcf,'B1_length_fig1.fig')

figure(2)
saveas(gcf,'B1_length_fig2.fig')