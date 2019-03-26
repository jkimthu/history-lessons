%% G. birth size vs time


%  Goal: is the stabilization timescale of V_birth different from that of
%        growth rate?

%        plot (1) birth volume vs time of birth
%        plot (2) division volume vs time of division
%        plot (3) added volume vs time of birth
%        plot (4) added volume vs time of division
%        plot (5) interdivision time vs time of birth
%        plot (6) interdivision time vs time of division


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes and added mass, using isDrop
%       d) calculate mean and plot condition data




%  Last edit: jen, 2019 March 25
%  Commit: add interdivision time as a variable to plot against time



%  OK let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = 13:15; % use corresponding dataIndex values


% 0. initialize parameters for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
binsPerHour = 30;


% 0. initialize data matrix for concatenating between experiments
cycleData = [];


%% Part 1. collect cell cycle data for each experiment in array


for e = 1:length(exptArray)
    
    
    % 1. collect experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    

    
    % 2. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    if strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
        % single upshift and downshift data only uses larger width thresh
    end
    load(filename,'D5','T');
    
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end expType filename
    
    
    
    % 4. for each condition of interest...
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        
        
        
        % 6. trim data to full cell cycles ONLY
        curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
        fullData = conditionData(curveID > 0,:);
        clear curveID
        
        
        
        % 7. isolate time, curve ids, volume, and birth event data (drop)
        curveFinder = getGrowthParameter(fullData,'curveFinder');      % curve Finder
        isDrop = getGrowthParameter(fullData,'isDrop');                % isDrop == 1 marks a birth event
        volume = getGrowthParameter(fullData,'volume');                % volume (Va)
        timestamps = getGrowthParameter(fullData,'timestamp')./3600;   % raw time (sec)
        clear fullData
        
        
        
        % 8. identify unique cell cycles by ID number
        cellCycles = curveFinder(isDrop == 1);
        
        
        
        % 9. from each unique cell cycle, collect data
        ccData = nan(length(cellCycles),7);
        for cc = 1:length(cellCycles)
            
            currentTimes = timestamps(curveFinder == cellCycles(cc));
            ccData(cc,1) = currentTimes(1);        % birth timestamp
            ccData(cc,2) = currentTimes(end);      % division timestamp
            ccData(cc,3) = (currentTimes(end)-currentTimes(1))*60; % interdivision time in min
            
            currentVolumes = volume(curveFinder == cellCycles(cc));
            ccData(cc,4) = currentVolumes(1);      % birth size
            ccData(cc,5) = currentVolumes(end);    % division size
            ccData(cc,6) = currentVolumes(end) - currentVolumes(1); % added size
            clear currentTimes currentVolumes
            
            ccData(cc,7) = condition;

        end
        clear cc volume timestamps isDrop curveFinder
        
        
        
        % 10. remove cell cycle data with interdivision times < 10 min
        tau = ccData(:,3);
        ccData_tau = ccData(tau > 10,:);
        clear tau
        
        
        
        % 11. remove cell cycle data with negative added volumes
        addedVol = ccData_tau(:,6);
        traits = ccData_tau(addedVol > 0,:);
        clear addedVol ccData_tau ccData
        
        
        
        % 12. remove data from post-bubble timestamps (no bubbles)
        maxTime = bubbletime(condition);
        birthTimes = traits(:,1);
        
        if maxTime > 0
            traits_final = traits(birthTimes <= maxTime,:);
        else
            traits_final = traits;
        end
        clear birthTimes maxTime cellCycles traits
        
        
        
        % 13. concatenate condition data with total cell cycle data
        cycleData = [cycleData; traits_final];
        size(traits_final)
        clear traits_final
        
    end
    
    
end
clear condition conditionData ans bubbletime date


%% Part 2. bin and plot data

for c = 1:4 % condition
        
        % 1. define plotting color
        color = rgb(palette(c));
        
        
        % 2. isolate condition data
        condVal = cycleData(:,7);
        cycleData_cond = cycleData(condVal == c,:);
        
        
        % 3. bin data into time bins
        timestamps_birth = cycleData_cond(:,1);
        timestamps_division = cycleData_cond(:,2);
        
        bins_div = ceil(timestamps_division * binsPerHour);
        bins_birth = ceil(timestamps_birth * binsPerHour);
        clear timestamps_birth timestamps_division
        
        
        % 1) birth volume vs time of birth
        volume_birth = cycleData_cond(:,4);
        volBirth_binnedByBT = accumarray(bins_birth,volume_birth,[],@(x) {x});
        
        % 2) division volume vs time of division
        volume_division = cycleData_cond(:,5);
        volDiv_binnedByDT = accumarray(bins_div,volume_division,[],@(x) {x});

        % 3) added volume vs time of birth
        volume_added = cycleData_cond(:,6);
        volAdded_binnedByBT = accumarray(bins_birth,volume_added,[],@(x) {x});
        
        % 4) added volume vs time of division
        volAdded_binnedByDT = accumarray(bins_div,volume_added,[],@(x) {x});
        clear volume_birth volume_division volume_added
        
        % 5) interdivision time vs time of birth
        interdivT = cycleData_cond(:,3);
        interdiv_binnedByBT = accumarray(bins_birth,interdivT,[],@(x) {x});
        
        % 6) interdivision time vs time of division
        interdiv_binnedByDT = accumarray(bins_div,interdivT,[],@(x) {x});
        clear interdivT
        
        
        % 14. generate time vectors for plotting
        timeVector_BT = (1:length(volBirth_binnedByBT))/binsPerHour;
        timeVector_DT = (1:length(volDiv_binnedByDT))/binsPerHour;
        
        
        
        
        
        % 15. plot

        figure(1) % birth volume vs. birth time
        errorbar(timeVector_BT,cellfun(@mean,volBirth_binnedByBT),cellfun(@std,volBirth_binnedByBT),'Color',color)
        hold on
        xlabel('time at birth (h)')
        ylabel('birth volume (cubic um)')
        title('stabilization of birth volume')
        axis([0 9.1 0 8])
        
        
        
        figure(2) % division volume vs. time at division
        errorbar(timeVector_DT,cellfun(@mean,volDiv_binnedByDT),cellfun(@std,volDiv_binnedByDT),'Color',color)
        hold on
        xlabel('time at division (h)')
        ylabel('division volume (cubic um)')
        title('stabilization of division volume')
        axis([0 9.1 0 16])
        
        
        
        figure(3) % added volume vs. time at birth
        errorbar(timeVector_BT,cellfun(@mean,volAdded_binnedByBT),cellfun(@std,volAdded_binnedByBT),'Color',color)
        hold on
        xlabel('time at birth (h)')
        ylabel('added volume (cubic um)')
        title('stabilization of added volume')
        axis([0 9.1 0 8])
        
        
        figure(4) % added volume vs. time at division
        errorbar(timeVector_DT,cellfun(@mean,volAdded_binnedByDT),cellfun(@std,volAdded_binnedByDT),'Color',color)
        hold on
        xlabel('time at division (h)')
        ylabel('added volume (cubic um)')
        title('stabilization of added volume')
        axis([0 9.1 0 8])
        
        
        figure(5) % interdiv time vs. time at birth
        errorbar(timeVector_BT,cellfun(@mean,interdiv_binnedByBT),cellfun(@std,interdiv_binnedByBT),'Color',color)
        hold on
        xlabel('time at birth (h)')
        ylabel('interdivision time (min)')
        title('stabilization of tau')
        axis([0 9.1 0 90])
        
        figure(6) % interdiv time vs. time at division
        errorbar(timeVector_DT,cellfun(@mean,interdiv_binnedByDT),cellfun(@std,interdiv_binnedByDT),'Color',color)
        hold on
        xlabel('time at division (h)')
        ylabel('interdivision time (min)')
        title('stabilization of tau')
        axis([0 9.1 0 90])
        
        clear timestamps_birth timestamps_division currentTimes cycleData_cond
       
        
end


