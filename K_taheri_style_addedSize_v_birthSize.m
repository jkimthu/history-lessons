%% Figure K. Taheri-style added size vs birth size

% goal: plot added size over binned birth size,
%       mirroring Fig 2a of Taheri-Araghi et al. Current Biology (2014)
%

% strategy:
%
%       0. initialize data
%       1. for all conditions in all fluctuating experiments
%       2. gather all added volumes and birth volumes for tracks...
%               - with birth after 3rd hour
%               - with division before bubble appearance (if any)
%       3. bin added volumes by birth volume
%       4. calculate mean and sem of added volumes
%       5. plot over birth volume bin


% versions of this plot:
%       v1. all cell cycles born after 3 hrs, dividing before any bubbles
%       v2. same as v1, but with x axis cropped to 0-10 cubic um
%       v3. all cycle cycles with negative added volumes removed
%       v4. all cell cycles born after FOUR hrs
%       v5. trim births before 3 hrs, normalize all values by condition mean
%           birth size
%       v6. with and without normalization, remove cell cycles with
%           V_division greater than 3x V_birth



% last update: jen, 2020 Jan 28
% commit: plot raw width instead of added width


% OK let's go!

%% initialize analysis

clc
clear

% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
experimentArray = [2:4,5:7,9:12,13:15];
experimentCount = length(experimentArray);


% 0. initialize environmental conditions
environment_order = {'low','3600','900','300','30','ave','high'};
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_ticks = zeros(length(environment_order),1);


% 0. initialize data collection
V_birth = cell(1,length(environment_order));
V_added = cell(1,length(environment_order));

L_birth = cell(1,length(environment_order));
L_added = cell(1,length(environment_order));

W_birth = cell(1,length(environment_order));
W_added = cell(1,length(environment_order));

T_birth = cell(1,length(environment_order));
T_division = cell(1,length(environment_order));

%% Part 1. loop through all experiments and compile data

for e = 1:experimentCount
    
    
    % 2. collect experiment meta data
    index = experimentArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    
    % 4. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    for condition = 1:length(bubbletime)
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);    % col 21 = condition vals
        
        
        
        % 6. trim data to full cell cycles ONLY
        curveID = getGrowthParameter(conditionData,'curveFinder');       % col 5  = curve finder (ID of curve in condition)
        fullData = conditionData(curveID > 0,:);
        clear curveID
        
        
        
        % 7. isolate time, curve ids, size (length, width, volume) and birth event data (drop)
        curveFinder = getGrowthParameter(fullData,'curveFinder');      % curve Finder
        isDrop = getGrowthParameter(fullData,'isDrop');                % isDrop == 1 marks a birth event
        lengthVals = getGrowthParameter(fullData,'length');            % length (um)
        width = getGrowthParameter(fullData,'width');                  % width
        volume = getGrowthParameter(fullData,'volume');                % volume (Va)
        timestamps = getGrowthParameter(fullData,'timestamp')./3600;   % raw time (sec)
        clear fullData
        
        

        
        % 8. prepare to assign condition data into a cell, where:
        %    column = environmental condition
        %    row = biological replicate
        
        % i. determine column no. of environmental condition
        if condition == 2
            eColumn = find(strcmp(environment_order,'low'));
        elseif condition == 3
            eColumn = find(strcmp(environment_order,'ave'));
        elseif condition == 4
            eColumn = find(strcmp(environment_order,'high'));
        else
            eColumn = find(strcmp(environment_order,num2str(timescale)));
        end
        environment_ticks(eColumn) = environment_ticks(eColumn) + 1;
        
        % ii. determine replicate no. of current condition data
        eRow = environment_ticks(eColumn);
        
        
        
        
        % 9. identify unique cell cycles by ID number
        cellCycles = curveFinder(isDrop == 1);
        
        
        
        % 10. from each unique cell cycle, collect data
        ccData = nan(length(cellCycles),12);
        for cc = 1:length(cellCycles)
            
            currentTimes = timestamps(curveFinder == cellCycles(cc));
            ccData(cc,1) = currentTimes(1);        % birth timestamp
            ccData(cc,2) = currentTimes(end);      % division timestamp
            ccData(cc,3) = (currentTimes(end)-currentTimes(1))*60; % interdivision time in min
            
            
            currentLengths = lengthVals(curveFinder == cellCycles(cc));
            ccData(cc,4) = currentLengths(1);      % birth length
            ccData(cc,5) = currentLengths(end);    % division length
            ccData(cc,6) = currentLengths(end) - currentLengths(1); % added length
            
            
            currentWidths = width(curveFinder == cellCycles(cc));
            ccData(cc,7) = currentWidths(1);      % birth width
            ccData(cc,8) = currentWidths(end);    % division width
            ccData(cc,9) = currentWidths(end) - currentWidths(1); % added width
            
            
            currentVolumes = volume(curveFinder == cellCycles(cc));
            ccData(cc,10) = currentVolumes(1);      % birth volume
            ccData(cc,11) = currentVolumes(end);    % division volume
            ccData(cc,12) = currentVolumes(end) - currentVolumes(1); % added volume
            clear currentTimes currentLengths currentWidths currentVolumes
            
            %ccData(cc,13) = condition;

        end
        clear cc lengthVals width volume timestamps isDrop curveFinder
        
        
        
        % 11. remove cell cycle data with interdivision times < 10 min
        tau = ccData(:,3);
        ccData_tau = ccData(tau > 10,:);
        clear tau
        
        
        
        % 12. remove cell cycle data with negative added lengths
        addedLength = ccData_tau(:,6);
        traits = ccData_tau(addedLength > 0,:);
        clear addedLength ccData_tau ccData
        
        
        
        % 13. remove data from post-bubble timestamps (no bubbles)
        maxTime = bubbletime(condition);
        birthTimes = traits(:,1);
        
        if maxTime > 0
            traits_final = traits(birthTimes <= maxTime,:);
        else
            traits_final = traits;
        end
        clear birthTimes maxTime cellCycles traits
        
        
        
        
        % 14. store condition data with total cell cycle data

        W_birth{eRow,eColumn} = traits_final(:,7);
        W_added{eRow,eColumn} = traits_final(:,9);
        
        L_birth{eRow,eColumn} = traits_final(:,4);
        L_added{eRow,eColumn} = traits_final(:,6);
        
        V_birth{eRow,eColumn} = traits_final(:,10);
        V_added{eRow,eColumn} = traits_final(:,12);
        
        T_birth{eRow,eColumn} = traits_final(:,1);
        T_division{eRow,eColumn} = traits_final(:,3);
        
    end
end
%%
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('addedSizeData.mat','V_birth','V_added','L_birth','L_added','W_birth','W_added','T_birth','T_division','environment_order','palette')

%% Part 2 - concatenate and plot: added volume/length vs birth volume/length


% 0. initialize data
clc
clear
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData.mat')



% 1. for each environment (column), concatenate data across replicates
for eCol = 1:length(environment_order)
    
    
    % 2. isolate column and initialize data
    hasData = cellfun(@isempty,V_birth(:,eCol));
    
    % volume
    currCol_birthvols = V_birth(hasData==0,eCol);
    currCol_addedvols = V_added(hasData==0,eCol);
    concat_Vbirth = [];
    concat_Vadded = [];
    
    % length
    currCol_birthlengths = L_birth(hasData==0,eCol);
    currCol_addedlengths = L_added(hasData==0,eCol);
    concat_Lbirth = [];
    concat_Ladded = [];
    
    % width
    currCol_birthwidths = W_birth(hasData==0,eCol);
    currCol_addedwidths = W_added(hasData==0,eCol);
    concat_Wbirth = [];
    concat_Wadded = [];
    clear hasData
    
    
    % 3. concatenate data by looping through replicates
    for replicate = 1:length(currCol_birthvols)
        
        % volume
        rep_Vbirths = currCol_birthvols{replicate};
        rep_Vadded = currCol_addedvols{replicate};
        concat_Vbirth = [concat_Vbirth; rep_Vbirths];
        concat_Vadded = [concat_Vadded; rep_Vadded];
        
        % length
        rep_Lbirths = currCol_birthlengths{replicate};
        rep_Ladded = currCol_addedlengths{replicate};
        concat_Lbirth = [concat_Lbirth; rep_Lbirths];
        concat_Ladded = [concat_Ladded; rep_Ladded];
        
        % width
        rep_Wbirths = currCol_birthwidths{replicate};
        rep_Wadded = currCol_addedwidths{replicate};
        concat_Wbirth = [concat_Wbirth; rep_Wbirths];
        concat_Wadded = [concat_Wadded; rep_Wadded];
        
        
    end
    
    
    
    % 5. remove unphysical events, where added volume is negative
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Lbirth = concat_Lbirth(concat_Vadded > 0);
    conceivables_Ladded = concat_Ladded(concat_Vadded > 0);
    
    conceivables_Wbirth = concat_Wbirth(concat_Vadded > 0);
    conceivables_Wadded = concat_Wadded(concat_Vadded > 0);
    clear concat_Vbirth concat_Vadded concat_Lbirth concat_Ladded concat_Wbirth concat_Wadded
    
    
    
    % 6. normalize all values by mean of birth size
    median_Vbirth = median(conceivables_Vbirth);
    median_Lbirth = median(conceivables_Lbirth);
    median_Wbirth = median(conceivables_Wbirth);
    
    
    
    % 7. bin birth sizes by every 0.1 cubic um 
    
    % volume
    birthBins_vol = floor(conceivables_Vbirth*10);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ 10;
    
    
    % length
    birthBins_length = floor(conceivables_Lbirth*10);
    bins = 1:max(birthBins_length);
    xbins_length = bins'/10;

    
    % width
    birthBins_width = floor(conceivables_Wbirth*10);
    bins = 1:max(birthBins_width);
    xbins_width = bins'/10;
    clear  bin 
    
    
    
    
    % 7. accumulate added size by birth bin
    
    % volume
    binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);

    
    % length
    binned_addedLengths = accumarray(birthBins_length,conceivables_Ladded,[],@(x) {x});
    binned_l_means = cellfun(@mean,binned_addedLengths);
    binned_l_stds = cellfun(@std,binned_addedLengths);
    binned_l_counts = cellfun(@length,binned_addedLengths);
    binned_l_sems = binned_l_stds./sqrt(binned_l_counts);
    

    % width
    binned_addedWidths = accumarray(birthBins_width,conceivables_Wadded,[],@(x) {x});
    binned_w_means = cellfun(@mean,binned_addedWidths);
    binned_w_stds = cellfun(@std,binned_addedWidths);
    binned_w_counts = cellfun(@length,binned_addedWidths);
    binned_w_sems = binned_w_stds./sqrt(binned_w_counts);
    
    
    
    
    % 8. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};
    
    
    % volume, raw
    figure(1)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_v_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_v_means,binned_v_sems,'.','Color',color)
    axis([0 10 0 10])
    title(condition)
    if eCol == 1
        ylabel('added volume,mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
   
    
    % length, raw
    figure(2)
    subplot(1,length(environment_order),eCol)
    plot(xbins_length,binned_l_means,'o','Color',color)
    hold on
    errorbar(xbins_length,binned_l_means,binned_l_sems,'.','Color',color)
    axis([0 6 0 6])
    title(condition)
    if eCol == 1
        ylabel('added length, mean + sem')
    end
    if eCol == 4
        xlabel('birth length, mean + sem')
    end
    
    
    
    % width, raw
    figure(3)
    subplot(1,length(environment_order),eCol)
    plot(xbins_width,binned_w_means,'o','Color',color)
    hold on
    errorbar(xbins_width,binned_w_means,binned_w_sems,'.','Color',color)
    axis([0.7 1.6 -0.2 0.2])
    title(condition)
    if eCol == 1
        ylabel('added length, mean + sem')
    end
    if eCol == 4
        xlabel('birth length, mean + sem')
    end
    
    
    
end


%% Part 3 - concatenate and plot: added length/width vs birth volume


% 0. initialize data
clc
clear
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData.mat')



% 1. for each environment (column), concatenate data across replicates
for eCol = 1:length(environment_order)
    
    
    % 2. isolate column and initialize data
    hasData = cellfun(@isempty,V_birth(:,eCol));
    
    % volume
    currCol_birthvols = V_birth(hasData==0,eCol);
    currCol_addedvols = V_added(hasData==0,eCol);
    concat_Vbirth = [];
    concat_Vadded = [];
    
    % length
    currCol_birthlengths = L_birth(hasData==0,eCol);
    currCol_addedlengths = L_added(hasData==0,eCol);
    concat_Lbirth = [];
    concat_Ladded = [];
    
    % width
    currCol_birthwidths = W_birth(hasData==0,eCol);
    currCol_addedwidths = W_added(hasData==0,eCol);
    concat_Wbirth = [];
    concat_Wadded = [];
    clear hasData
    
    
    % 3. concatenate data by looping through replicates
    for replicate = 1:length(currCol_birthvols)
        
        % volume
        rep_Vbirths = currCol_birthvols{replicate};
        rep_Vadded = currCol_addedvols{replicate};
        concat_Vbirth = [concat_Vbirth; rep_Vbirths];
        concat_Vadded = [concat_Vadded; rep_Vadded];
        
        % length
        rep_Lbirths = currCol_birthlengths{replicate};
        rep_Ladded = currCol_addedlengths{replicate};
        concat_Lbirth = [concat_Lbirth; rep_Lbirths];
        concat_Ladded = [concat_Ladded; rep_Ladded];
        
        % width
        rep_Wbirths = currCol_birthwidths{replicate};
        rep_Wadded = currCol_addedwidths{replicate};
        concat_Wbirth = [concat_Wbirth; rep_Wbirths];
        concat_Wadded = [concat_Wadded; rep_Wadded];
        
        
    end
    
    
    
    % 5. remove unphysical events, where added volume is negative
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Lbirth = concat_Lbirth(concat_Vadded > 0);
    conceivables_Ladded = concat_Ladded(concat_Vadded > 0);
    
    conceivables_Wbirth = concat_Wbirth(concat_Vadded > 0);
    conceivables_Wadded = concat_Wadded(concat_Vadded > 0);
    clear concat_Vbirth concat_Vadded concat_Lbirth concat_Ladded concat_Wbirth concat_Wadded
    
    
    
    
    % 6. normalize all values by mean of birth size
    median_Vbirth = median(conceivables_Vbirth);

    
    
    
    % 7. bin birth sizes by every 0.1 cubic um 
    
    % volume
    birthBins_vol = floor(conceivables_Vbirth*10);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ 10;
    clear  bin 
    
    
    
    
    % 8. accumulate added size by birth bin
    
    % volume
    binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);

    
    % length
    binned_addedLengths = accumarray(birthBins_vol,conceivables_Ladded,[],@(x) {x});
    binned_l_means = cellfun(@mean,binned_addedLengths);
    binned_l_stds = cellfun(@std,binned_addedLengths);
    binned_l_counts = cellfun(@length,binned_addedLengths);
    binned_l_sems = binned_l_stds./sqrt(binned_l_counts);
    

    % width
    binned_addedWidths = accumarray(birthBins_vol,conceivables_Wadded,[],@(x) {x});
    binned_w_means = cellfun(@mean,binned_addedWidths);
    binned_w_stds = cellfun(@std,binned_addedWidths);
    binned_w_counts = cellfun(@length,binned_addedWidths);
    binned_w_sems = binned_w_stds./sqrt(binned_w_counts);
    
    
    
    
    % 9. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};
    
    
    % volume, raw
    figure(1)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_v_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_v_means,binned_v_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,0,'o','Color',color)
    axis([0 10 0 10])
    title(condition)
    if eCol == 1
        ylabel('added volume,mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
   
    
    % length, raw
    figure(2)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_l_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_l_means,binned_l_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,0,'o','Color',color)
    axis([0 10 0 6])
    title(condition)
    if eCol == 1
        ylabel('added length, mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    
    
    % width, raw
    figure(3)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_w_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_w_means,binned_w_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,-0.5,'o','Color',color)
    axis([0 10 -0.5 0.5])
    title(condition)
    if eCol == 1
        ylabel('added width, mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    
    % pdf
    count_sum = sum(binned_w_counts);
    densityFunc = binned_w_counts./count_sum;
    
    figure(4)
    subplot(1,length(environment_order),eCol)
    bar(xbins_vol,densityFunc,'FaceColor',color)
    hold on
    plot(median_Vbirth,0,'o','Color',rgb('SlateGray'))
    hold on
    txt = strcat('n=',num2str(count_sum));
    text(2,0.13,txt)
    axis([0 10 0 0.12])
    title(condition)
    if eCol == 1
        ylabel('pdf')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    
end

% 10. save plots
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/K_addedSize_v_birthSize')

figure(1)
plotName = strcat('K-addedVol_v_birthVol');
saveas(gcf,plotName,'epsc')

figure(2)
plotName = strcat('K-addedLength_v_birthVol');
saveas(gcf,plotName,'epsc')

figure(3)
plotName = strcat('K-addedWidth_v_birthVol');
saveas(gcf,plotName,'epsc')

figure(4)
plotName = strcat('K-pdf_v_birthVol');
saveas(gcf,plotName,'epsc')


%% Part 4 - concatenate and plot: added v, added l, birth width vs birth volume


% 0. initialize data
clc
clear
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData.mat')



% 1. for each environment (column), concatenate data across replicates
for eCol = 1:length(environment_order)
    
    
    % 2. isolate column and initialize data
    hasData = cellfun(@isempty,V_birth(:,eCol));
    
    % volume
    currCol_birthvols = V_birth(hasData==0,eCol);
    currCol_addedvols = V_added(hasData==0,eCol);
    concat_Vbirth = [];
    concat_Vadded = [];
    
    % length
    currCol_birthlengths = L_birth(hasData==0,eCol);
    currCol_addedlengths = L_added(hasData==0,eCol);
    concat_Lbirth = [];
    concat_Ladded = [];
    
    % width
    currCol_birthwidths = W_birth(hasData==0,eCol);
    %currCol_addedwidths = W_added(hasData==0,eCol);
    concat_Wbirth = [];
    concat_Wadded = [];
    clear hasData
    
    
    % 3. concatenate data by looping through replicates
    for replicate = 1:length(currCol_birthvols)
        
        % volume
        rep_Vbirths = currCol_birthvols{replicate};
        rep_Vadded = currCol_addedvols{replicate};
        concat_Vbirth = [concat_Vbirth; rep_Vbirths];
        concat_Vadded = [concat_Vadded; rep_Vadded];
        
        % length
        rep_Lbirths = currCol_birthlengths{replicate};
        rep_Ladded = currCol_addedlengths{replicate};
        concat_Lbirth = [concat_Lbirth; rep_Lbirths];
        concat_Ladded = [concat_Ladded; rep_Ladded];
        
        % width
        rep_Wbirths = currCol_birthwidths{replicate};
        %rep_Wadded = currCol_addedwidths{replicate};
        concat_Wbirth = [concat_Wbirth; rep_Wbirths];
        %concat_Wadded = [concat_Wadded; rep_Wadded];
        
        
    end
    
    
    
    % 5. remove unphysical events, where added volume is negative
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Lbirth = concat_Lbirth(concat_Vadded > 0);
    conceivables_Ladded = concat_Ladded(concat_Vadded > 0);
    
    conceivables_Wbirth = concat_Wbirth(concat_Vadded > 0);
    %conceivables_Wadded = concat_Wadded(concat_Vadded > 0);
    clear concat_Vbirth concat_Vadded concat_Lbirth concat_Ladded concat_Wbirth concat_Wadded
    
    
    
    
    % 6. normalize all values by mean of birth size
    median_Vbirth = median(conceivables_Vbirth);

    
    
    
    % 7. bin birth sizes by every 0.1 cubic um 
    
    % volume
    birthBins_vol = floor(conceivables_Vbirth*10);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ 10;
    clear  bin 
    
    
    
    
    % 8. accumulate size data by birth bin
    
    % added volume
    binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);

    
    % added length
    binned_addedLengths = accumarray(birthBins_vol,conceivables_Ladded,[],@(x) {x});
    binned_l_means = cellfun(@mean,binned_addedLengths);
    binned_l_stds = cellfun(@std,binned_addedLengths);
    binned_l_counts = cellfun(@length,binned_addedLengths);
    binned_l_sems = binned_l_stds./sqrt(binned_l_counts);
    

    % birth width
    binned_bWidths = accumarray(birthBins_vol,conceivables_Wbirth,[],@(x) {x});
    binned_w_means = cellfun(@mean,binned_bWidths);
    binned_w_stds = cellfun(@std,binned_bWidths);
    binned_w_counts = cellfun(@length,binned_bWidths);
    binned_w_sems = binned_w_stds./sqrt(binned_w_counts);
    
    
    
    
    % 9. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};
    
    
    % volume, raw
    figure(1)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_v_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_v_means,binned_v_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,0,'o','Color',color)
    axis([0 10 0 10])
    title(condition)
    if eCol == 1
        ylabel('added volume,mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
   
    
    % length, raw
    figure(2)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_l_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_l_means,binned_l_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,0,'o','Color',color)
    axis([0 10 0 6])
    title(condition)
    if eCol == 1
        ylabel('added length, mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    
    
    % width, raw
    figure(3)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_w_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_w_means,binned_w_sems,'.','Color',color)
    hold on
    plot(median_Vbirth,-0.5,'o','Color',color)
    axis([0 10 0.8 1.8])
    title(condition)
    if eCol == 1
        ylabel('birth width, mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    
    % pdf
%     count_sum = sum(binned_w_counts);
%     densityFunc = binned_w_counts./count_sum;
%     
%     figure(4)
%     subplot(1,length(environment_order),eCol)
%     bar(xbins_vol,densityFunc,'FaceColor',color)
%     hold on
%     plot(median_Vbirth,0,'o','Color',rgb('SlateGray'))
%     hold on
%     txt = strcat('n=',num2str(count_sum));
%     text(2,0.13,txt)
%     axis([0 10 0 0.12])
%     title(condition)
%     if eCol == 1
%         ylabel('pdf')
%     end
%     if eCol == 4
%         xlabel('birth volume, mean + sem')
%     end
%     
    
end

% 10. save plots
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/K_addedSize_v_birthSize')

figure(1)
plotName = strcat('K-addedVol_v_birthVol');
%saveas(gcf,plotName,'epsc')

figure(2)
plotName = strcat('K-addedLength_v_birthVol');
%saveas(gcf,plotName,'epsc')

figure(3)
plotName = strcat('K-addedWidth_v_birthVol');
%saveas(gcf,plotName,'epsc')

figure(4)
plotName = strcat('K-pdf_v_birthVol');
%saveas(gcf,plotName,'epsc')

