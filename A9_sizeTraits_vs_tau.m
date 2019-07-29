%% Figure A9: size traits vs interdivision time, heated scatter and bubble


%  Goal: split cell cycles by signal experienced, plot mean with standard
%        deviation in both directions



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell interdivision time and instantaneous growth rates
%  Part 3. visualize
%  Part 4. fit best line



%  Last edit: jen, 2019 July 27
%  Commit: add visualizations


%  OK let's go!



%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. create array of experiments to use in analysis
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
compiled_tau = cell(length(exptArray),1);
compiled_addedSize = cell(length(exptArray),1);
compiled_divSize = cell(length(exptArray),1);
compiled_signal = cell(length(exptArray),1);


%% Part 2. collect single cell interdivision times and instantaneous growth rates

%  Strategy:
%
%       1.  loop through each experiment to collect data
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. for each condition in experiment...
%                       5. isolate condition specific data
%                       6. to calculate growth rate, isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. trim condition and growth rate data to include only full cell cycles
%                       9. isolate isDrop, timestamp, and timeSinceBirth data
%                      10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
%                      11. remove zeros, which occur if no full track data exists at a drop
%                      12. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      13. truncate data to stabilized regions
%                      14. if no div data in steady-state, skip condition
%                          else, trim outliers (those 3 std dev away from median) from final dataset
%                      15. bin growth rates and signals by cell cycle, to match organization of birth size data
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
    interDiv = cell(length(bubbletime),1);
    delta = cell(length(bubbletime),1);
    Sdiv = cell(length(bubbletime),1);
    nSignal = cell(length(bubbletime),1);

    
    
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
        
        
        
        % 6. to calculate growth rate, isolate volume (Va), timestamp, drop, curve, and trackNum data
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
        
        
        
        % 9. calculate binary nutrient signals
        [binaryNutrientSignal, nScore] = nutrientScore(timescale,conditionData_fullOnly);
        
        
        
        % 10. for cell cycle measurements, isolate timestamp, isDrop, size and interdiv time data 
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        size = getGrowthParameter(conditionData_fullOnly,'volume');        % calculated va_vals (cubic um)
        addedSize = getGrowthParameter(conditionData_fullOnly,'addedVA');
        tau = getGrowthParameter(conditionData_fullOnly,'curveDurations'); % tau repeated for entire growth curve
        clear timestamps
        
        
        
        % 11. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        isBirth = find(isDrop == 1);
        convert = isBirth - 1;
        isDiv = [convert(2:end); length(size)];
        
        final_birthSize = size(isDrop==1);
        final_dsize = size(isDiv);
        final_addedSize = addedSize(isDrop==1);
        final_tau = tau(isDrop==1);
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        clear conditionData_fullOnly timestamps_hr convert isDrop isDiv isBirth
        
        

        % 12. remove zeros, which occur if no full track data exists at a drop
        dSize = final_dsize(final_birthSize > 0);
        aSize = final_addedSize(final_birthSize > 0);
        interdivisionT = final_tau(final_birthSize > 0);
        birthTimestamps = finalTimestamps(final_birthSize > 0);
        curveIDs = curveIDs_unique(final_birthSize > 0);
        clear final_dsize final_addedSize finalTimestamps addedSize size final_tau tau
        
        
        
        % 13. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            interdivisionT_bubbleTrimmed = interdivisionT(birthTimestamps <= maxTime,:);
            dSize_bubbleTrimmed = dSize(birthTimestamps <= maxTime,:);
            aSize_bubbleTrimmed = aSize(birthTimestamps <= maxTime,:);
            curveIDs_bubbleTrimmed_cc = curveIDs(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
             
        else
            
            interdivisionT_bubbleTrimmed = interdivisionT;
            dSize_bubbleTrimmed = dSize;
            aSize_bubbleTrimmed = aSize;
            curveIDs_bubbleTrimmed_cc = curveIDs;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear maxTime interdivisionT dSize aSize curveIDs birthTimestamps
        
        
        
        % 14. truncate data to stabilized regions
        minTime = 3;
        interdivisionT_fullyTrimmed = interdivisionT_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        dSize_fullyTrimmed = dSize_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        aSize_fullyTrimmed = aSize_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        curveIDs_fullyTrimmed_cc = curveIDs_bubbleTrimmed_cc(birthTimestamps_bubbleTrimmed >= minTime,:);        
        clear interdivisionT_bubbleTrimmed  dSize_bubbleTrimmed aSize_bubbleTrimmed curveIDs_bubbleTrimmed_cc birthTimestamps_bubbleTrimmed
        
        
        
        
        % 15. if no div data in steady-state, skip condition
        if isempty(interdivisionT_fullyTrimmed) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            tau_median = median(interdivisionT_fullyTrimmed);
            tau_std_temp = std(interdivisionT_fullyTrimmed);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            tau_temp = interdivisionT_fullyTrimmed(interdivisionT_fullyTrimmed <= (tau_median+tau_std_temp*3)); % cut largest vals, over 3 std out
            dSize_temp = dSize_fullyTrimmed(dSize_fullyTrimmed <= (tau_median+tau_std_temp*3));
            aSize_temp = aSize_fullyTrimmed(aSize_fullyTrimmed <= (tau_median+tau_std_temp*3));
            IDs_temp = curveIDs_fullyTrimmed_cc(interdivisionT_fullyTrimmed <= (tau_median+tau_std_temp*3));
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            tau_final = tau_temp(tau_temp >= (tau_median-tau_std_temp*3));          % cut smallest vals, over 3 std out
            dSize_final = dSize_temp(tau_temp >= (tau_median-tau_std_temp*3)); 
            aSize_final = aSize_temp(tau_temp >= (tau_median-tau_std_temp*3)); 
            IDs_final = IDs_temp(tau_temp >= (tau_median-tau_std_temp*3));   
            clear tau_median tau_std_temp tau_temp IDs_temp dSize_temp aSize_temp
            
            % iv. remove corresponding growth rates and signals from datasets
            trimmedIDs = setdiff(curveIDs_unique,IDs_final);    % curve IDs in growth rate dataset, NOT in final IDs trimmed by cell cycle
            toTrim = ismember(curveIDs_fullOnly,trimmedIDs);   % vector of what to trim or not in growth rate
            trimmed_curves_insta = curveIDs_fullOnly(toTrim == 0);
            trimmed_mus = growthRates_fullOnly(toTrim == 0);
            trimmed_signals = binaryNutrientSignal(toTrim == 0);
            clear toTrim trimmedIDs curveIDs_fullOnly growthRates_fullOnly binaryNutrientSignal
            
            
                 
            % 16. bin growth rates and nutrient signal by cell cycle, to match organization of birth size data
            mus_binned = accumarray(trimmed_curves_insta,trimmed_mus,[],@(x) {x});
            signals_binned = accumarray(trimmed_curves_insta,trimmed_signals,[],@(x) {x});
            mus = mus_binned(~cellfun('isempty',mus_binned));
            signals = signals_binned(~cellfun('isempty',signals_binned));
            clear trimmed_curves_insta trimmed_mus trimmed_signals
            
            
            
            % 16. remove cell cycles with negative mean growth rate
            meanGR = cellfun(@nanmean,mus);
            
            mus = mus(meanGR > 0);
            tau_final = tau_final(meanGR > 0);
            dSize_final = dSize_final(meanGR > 0);
            aSize_final = aSize_final(meanGR > 0);
            signals = signals(meanGR > 0);
            
            
            % 17. store condition data into one variable per experiment
            interDiv{condition} = tau_final;
            delta{condition} = aSize_final;
            Sdiv{condition} = dSize_final;
            nSignal{condition} = signals;
        
            
        end
        
    end
      
    
    % 18. store experiment data into single variable for further analysis
    compiled_tau{e} = interDiv;
    compiled_addedSize{e} = delta;
    compiled_divSize{e} = Sdiv;
    compiled_signal{e} = nSignal;
    
end

% 19. save hard earned data
save('A9_data_volume.mat','compiled_tau','compiled_addedSize','compiled_divSize','compiled_signal','exptArray')

%% Part 3. data visualization


clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A9_data_volume.mat')
load('A6_data.mat','compiled_lambda')
condition = 1;


% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','ForestGreen','Orange','GoldenRod'};
shape = 'o';


% 0 . initialize plot limits
plot_tau = 90;
plot_delta = 6;
plot_divSize = 11;
max_lambda = 5;


% 0. initialize signal classifications
onlyH = 1;
onlyL = 2;
H2L = 3;
L2H = 4;
HLH = 5;
LHL = 6;
HLHL =7;
LHLH = 8;
sigArray = [onlyH; onlyL; H2L; L2H; HLH; LHL; HLHL; LHLH];
sigNames = {'onlyH','onlyL','H2L','L2H','HLH','LHL','HLHL','LHLH'};




% 1. loop through experiments to format data and plot
for ee = 11:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(ee);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    clear index
    
    
    % 3. isolate experiment data
    eeTaus = compiled_tau{ee}{condition}./60; % convert from sec to h
    eeMus = compiled_lambda{ee}{condition};
    eeSignals = compiled_signal{ee}{condition};
    eeDeltas = compiled_addedSize{ee}{condition};
    eeDivSizes = compiled_divSize{ee}{condition};
    
    
    % 4. calculate lambda and nScores from compiled mus
    eeLambdas = cellfun(@nanmean,eeMus);
    clear eeMus
    
    
    % 5. remove data with too large growth rates
    taus = eeTaus(eeLambdas <= max_lambda);
    nSignals = eeSignals(eeLambdas <= max_lambda);
    lambdas = eeLambdas(eeLambdas <= max_lambda);
    deltas = eeDeltas(eeLambdas <= max_lambda);
    divSizes = eeDivSizes(eeLambdas <= max_lambda);
    clear eeNscores eeLambdas eeSignals eeTaus eeDeltas eeDivSizes
    
    
    % 6. for each curve, determine signal type
    signalType = zeros(length(nSignals),1);
    if timescale == 3600
        
        for cc = 1:length(nSignals)
            
            currSignal = nSignals{cc};
            ds_dt = diff(currSignal);
            currShifts = length(find(ds_dt ~= 0));
            
            if currShifts == 0
                if currSignal(1) == 1
                    signalType(cc) = onlyH;
                else
                    signalType(cc) = onlyL;
                end
            elseif currShifts == 1
                if currSignal(1) == 1
                    signalType(cc) = H2L;
                else
                    signalType(cc) = L2H;
                end
            elseif currShifts == 2
                if currSignal(1) == 1
                    signalType(cc) = HLH;
                else
                    signalType(cc) = LHL;
                end
            elseif currShifts == 3
                if currSignal(1) == 1
                    signalType(cc) = HLHL;
                else
                    signalType(cc) = LHLH;
                end
            end
        end
        
    else
        error('wrong signal classifier for current timescale!')
    end
    clear currSignal ds_dt currShifts cc
    
    
    
    % 7. accumulate cell cycle phenotypes by signal class
    tau_binned = accumarray(signalType,taus,[],@(x) {x});
    delta_binned = accumarray(signalType,deltas,[],@(x) {x});
    divSize_binned = accumarray(signalType,divSizes,[],@(x) {x});
    
    
    % 8. calculate mean and standard error for each signal class
    tau_mean = cellfun(@mean,tau_binned);
    tau_std = cellfun(@std,tau_binned);
    tau_count = cellfun(@length,tau_binned);
    tau_sem = tau_std./sqrt(tau_count);
    
    delta_mean = cellfun(@mean,delta_binned);
    delta_std = cellfun(@std,delta_binned);
    delta_count = cellfun(@length,delta_binned);
    delta_sem = delta_std./sqrt(delta_count);
    
    divSize_mean = cellfun(@mean,divSize_binned);
    divSize_std = cellfun(@std,divSize_binned);
    divSize_count = cellfun(@length,divSize_binned);
    divSize_sem = divSize_std./sqrt(divSize_count);
    
    
    % 9. plot mean and standard error of each signal class
    figure(1)
    for sc = 1:length(delta_mean)
        color = rgb(palette{sc});
        if isempty(delta_mean(sc)) == 1
            continue
        else
            errorbar(delta_mean(sc),tau_mean(sc),tau_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            hold on
            errorbar(delta_mean(sc),tau_mean(sc),delta_sem(sc),'horizontal','Color',color)
            
            %errorbar(tau_mean(sc),delta_mean(sc),delta_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            %hold on
            %errorbar(tau_mean(sc),delta_mean(sc),tau_sem(sc),'horizontal','Color',color)
        end
    end
    ylabel('Interdivison time')
    xlabel('Added size')
    %axis([0 plot_delta 0 plot_tau])
    
    figure(2)
    for sc = 1:length(delta_mean)
        color = rgb(palette{sc});
        if isempty(delta_mean(sc)) == 1
            continue
        else
            errorbar(divSize_mean(sc),tau_mean(sc),tau_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            hold on
            errorbar(divSize_mean(sc),tau_mean(sc),divSize_sem(sc),'horizontal','Color',color)
            
        end
    end
    ylabel('Interdivison time')
    xlabel('Division size')
    %axis([0 plot_divSize 0 plot_tau])
    clear color sc tau_mean tau_sem tau_std tau_count delta_mean delta_sem delta_std delta_count
    clear divSize_mean divSize_sem divSize_std divSize_count
    clear signalType taus nSignals nScores delta_binned delta tau_binned
    
end
clear onlyH onlyL H2L L2H HLH LHL HLHL LHLH


% 10. plot expected tau = 1/lambda relationship
% x = linspace(0.1,plot_delta,100);
% y = 60./x; % convert from h to min
% 
% figure(1)
% hold on
% plot(x,y,'Color',rgb('SlateGray'))


%%
% 11. save plot
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('A9-delta-v-tau-bubbled');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(2)
plotName = strcat('A9-divSize-v-tau-bubbled');
saveas(gcf,plotName,'epsc')
close(gcf)
