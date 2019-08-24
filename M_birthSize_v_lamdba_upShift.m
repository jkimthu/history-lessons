%% M. birth size vs lambda for increments after a single shift


%  Goal: what is the path between two steady-states?

%        plot birth volume vs mean growth rate from data just before
%        and for as long as possible after a single shift


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes and added mass, using isDrop
%       d) calculate mean and plot condition data




%  Last edit: jen, 2019 Aug 23
%  Commit: re-do such that upshift and downshift have same axis



%  OK let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = 21:22; % use corresponding dataIndex values


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. initialize parameters for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
shape = 'o';
binsPerHour = 2;


% 0. initialize data matrix for concatenating between experiments
cycleData = [];


%% Part 1. collect cell cycle data for each experiment in array


for e = 1:length(exptArray)
    
    
    % 1. collect experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    shiftTime = storedMetaData{index}.shiftTime/3600; % sec convert to h
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
    
    
    
    % 3. for each condition of interest...
    condArray = [1,2,3,4];
    for ii = 1:length(condArray)
        
        condition = condArray(ii);
        
        
        % 4. compile experiment data matrix
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end
        
        
        
        % 5. for growth rate calculations, isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
        
        
        
        % 6. calculate growth rate
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes isDrop trackNum timestamps_sec
        

        
        % 7. trim data to full cell cycles ONLY
        curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
        fullData = conditionData(curveID > 0,:);
        fullGR = growthRates(curveID > 0,:);
        clear curveID growthRates_all growthRates
        
        
        
        % 8. isolate time, curve ids, volume, and birth event data (drop)
        curveFinder = getGrowthParameter(fullData,'curveFinder');      % curve Finder
        isDrop = getGrowthParameter(fullData,'isDrop');                % isDrop == 1 marks a birth event
        volume = getGrowthParameter(fullData,'volume');                % volume (Va)
        timestamps = getGrowthParameter(fullData,'timestamp')./3600;   % raw time (sec) converted to h
        clear fullData
        
        
        
        % 9. identify unique cell cycles by ID number
        cellCycles = curveFinder(isDrop == 1);
        
        
        
        % 10. from each unique cell cycle, collect data
        ccData = nan(length(cellCycles),9);
        for cc = 1:length(cellCycles)
            
            currentTimes = timestamps(curveFinder == cellCycles(cc));
            ccData(cc,1) = currentTimes(1) - shiftTime;        % birth timestamp relative to shift
            ccData(cc,2) = currentTimes(end) - shiftTime;      % division timestamp relative to shift
            ccData(cc,3) = (currentTimes(end)-currentTimes(1))*60; % interdivision time in min
            
            currentVolumes = volume(curveFinder == cellCycles(cc));
            ccData(cc,4) = currentVolumes(1);      % birth size
            ccData(cc,5) = currentVolumes(end);    % division size
            ccData(cc,6) = currentVolumes(end) - currentVolumes(1); % added size
            clear currentTimes currentVolumes
            
            ccData(cc,7) = condition;
            ccData(cc,8) = e;
            
            currentMus = fullGR(curveFinder == cellCycles(cc));
            ccData(cc,9) = nanmean(currentMus);

        end
        clear cc volume timestamps isDrop curveFinder currentMus currentVolumes currentTimes
        
        
        
        % 11. remove cell cycle data with interdivision times < 9 min
        tau = ccData(:,3);
        ccData_tau = ccData(tau > 9,:);
        clear tau
        
        
        
        % 12. remove cell cycle data with negative added volumes
        addedVol = ccData_tau(:,6);
        traits = ccData_tau(addedVol > 0,:);
        clear addedVol ccData_tau ccData
        
        
        
        % 13. remove data from post-bubble timestamps (no bubbles)
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
clear D5 T ii index shiftTime xys timescale

% 14. save data
save('M_data_upshift.mat','cycleData','condArray','shape','palette','binsPerHour')


%% Part 2. bin and plot data
clear
clc
load('M_data_upshift.mat')

% time bins of interest in this analysis 
minBin = 20; % 18-20 are preshift bins
maxBin = 30; % final bin before weird end behavior starts

for jj = 1:length(condArray) % condition
        
    
        % 1. determine current condition of interest
        c = condArray(jj);
    
      
        % 2. isolate condition data
        condVal = cycleData(:,7);
        cycleData_cond = cycleData(condVal == c,:);
        
        
        % 3. bin data into time bins
        %    note: because timestamps are relative to shift time, the data
        %    must first be shifted to avoided negative values, and then
        %    shifted back after binning.
        
        timestamps_birth = cycleData_cond(:,1) + 10;
        timestamps_division = cycleData_cond(:,2) + 10;
        bins_div = ceil(timestamps_division * binsPerHour);
        bins_birth = ceil(timestamps_birth * binsPerHour);
        clear timestamps_birth timestamps_division
        
        
        % birth volume vs time of birth
        volume_birth = cycleData_cond(:,4);
        volBirth_binnedByBT = accumarray(bins_birth,volume_birth,[],@(x) {x});       
  
        % lambda vs time of birth
        lambda = cycleData_cond(:,9);
        lambda_bindedByBT = accumarray(bins_birth,lambda,[],@(x) {x});
        
        
        
        % 4. generate time vectors for plotting
        %    note: in time vector, zero = bin immediately before shift
        timeVector_BT = ((1:length(volBirth_binnedByBT))/binsPerHour) - 10;
        
        
        
        % 5. plot sanity check
%         figure(1) % birth volume vs. birth time
%         errorbar(timeVector_BT,cellfun(@mean,volBirth_binnedByBT),cellfun(@std,volBirth_binnedByBT),'Color',color)
%         hold on
%         xlabel('time at birth (h)')
%         ylabel('birth volume (cubic um)')
%         title('stabilization of birth volume')
        


        % 6. plot what we want for reals
        figure(1) % birth volume vs. lambda
        if c == 1
            
            % i. isolate data to time bins of interest
            lambda_final = lambda_bindedByBT(minBin:maxBin);
            volBirth_final = volBirth_binnedByBT(minBin:maxBin);
            timeVector_final = timeVector_BT(minBin:maxBin);
        
                    
            % ii. single shift in 30 min bins
            shift_Vb = cellfun(@mean,volBirth_final);
            shift_Vb_std = cellfun(@std,volBirth_final);
            shift_Vb_count = cellfun(@length,volBirth_final);
            shift_Vb_sem = shift_Vb_std./sqrt(shift_Vb_count);
            
            shift_gr = cellfun(@mean,lambda_final);
            shift_gr_std = cellfun(@std,lambda_final);
            shift_gr_count = cellfun(@length,lambda_final);
            shift_gr_sem = shift_gr_std./sqrt(shift_gr_count);
            
            cmap = parula(length(shift_gr)); % spread parula over length of 2d vector.
            for tt = 1:length(shift_gr)
                hold on
                plot(shift_gr(tt),log(shift_Vb(tt)),'Color',cmap(tt,:),'Marker',shape,'MarkerSize',10,'LineWidth',2)
                colorbar
            end
            clear shift_Vb shift_Vb_std shift_Vb_count shift_Vb_sem shift_gr shift_gr_std shift_gr_count shift_gr_sem
            clear tt
            
        else
            
            % i. define plotting color
            color = rgb(palette(c));
     

            % ii. isolate data from bins
            lambda_trim1 = lambda(bins_birth >= minBin);
            volume_birth_trim1 = volume_birth(bins_birth >= minBin);
            bins_birth_trim1 = bins_birth(bins_birth >= minBin);
            
            lambda_trim2 = lambda_trim1(bins_birth_trim1 <= maxBin);
            volume_birth_trim2 = volume_birth_trim1(bins_birth_trim1 <= maxBin);
            %bins_birth_trim2 = bins_birth_trim1(bins_birth_trim1 <= maxBin);
            
            
            % iii. calculate mean and sem
            lambda_steady(c) = mean(lambda_trim2);
            volume_birth_steady(c) = mean(volume_birth_trim2);
            
            % iv. plot
            figure(1)
            hold on
            plot(lambda_steady(c),log(volume_birth_steady(c)),'Color',color,'MarkerFaceColor',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
            
            
        end
        clear timestamps_birth timestamps_division currentTimes cycleData_cond    
        
end

% 7. plot fit line between steady environments
lambda_steady = lambda_steady(2:4);
volume_birth_steady = volume_birth_steady(2:4);

fit = polyfit(lambda_steady,log(volume_birth_steady),1);

x = linspace(lambda_steady(1),lambda_steady(end),10);
y = fit(1)*x + fit(2);

figure(1)
hold on
plot(x,y,'Color',rgb('SlateGray'))
axis([0.3 3.8 0.5 2])
xlabel('mean mu')
ylabel('ln(mean birth volume)')
title('single upshifts compiled')

%%
