%% Figure B3: tau vs lambda, heated scatter


%  Goal: plot a heated scatter of single-cell interdivision time vs mean growth rate.
%        consider only cell cycles born after 3 hrs


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, division volumes and time at birth
%       d) plot raw condition data and heatmap version, for each cycle
%       e) save both figures per experiment



%  Last edit: jen, 2019 July 17
%  Commit: edit to analyze all timescales on windows workstation



%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment = {'fluc','low','ave','high'};

%%
% 1. for all experiments in dataset
ts_all = {2:4;5:7;9:12;13:15};
%exptArray = 13:15;

for ts = 1:length(ts_all)
    
    % 0. initialize array of experiments to use in analysis, then loop through each
    exptArray = ts_all{ts}; % use corresponding dataIndex values
    
    for condition = 1:4 
        
        % 0. initialize array for concatenation of final cell cycles from each experimental dataset
        lambda = [];
        tau = [];
        cc_lengths = [];
        t_division = [];
        
        
        % 1. for all experiments in dataset
        for e = 1:length(exptArray)
            
            % 1. collect experiment meta data
            index = exptArray(e);
            date = storedMetaData{index}.date;
            timescale = storedMetaData{index}.timescale;
            expType = storedMetaData{index}.experimentType;
            bubbletime = storedMetaData{index}.bubbletime;
            xys = storedMetaData{index}.xys;
            disp(strcat(date, ': analyze!'))
            
            
            % 2. load measured data
            %experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
            %cd(experimentFolder)
            filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
            load(filename,'D5','T');
            
            
            % 3. compile experiment data matrix
            xy_start = xys(condition,1);
            xy_end = xys(condition,end);
            conditionData = buildDM(D5,T,xy_start,xy_end,index,expType);
            %exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
            clear D5 T xy_start xy_end e
            
            
            
            %for condition = 1:length(environment)
            
            
            % 5. isolate condition specific data
            %conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
            
            
            % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
            volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
            timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
            isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
            curveID = getGrowthParameter(conditionData,'curveFinder');       % curve finder (ID of curve in condition)
            trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
            
            
            
            % 7. calculate growth rate
            growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveID,trackNum);
            growthRates = growthRates_all(:,specificColumn);
            clear volumes timestamps_sec isDrop trackNum
            
            
            
            % 8. trim data to full cell cycles ONLY
            conditionData_fullOnly = conditionData(curveID > 0,:);
            growthRates_fullOnly = growthRates(curveID > 0,:);
            clear curveID conditionData growthRates growthRates_all
            
            
            
            % 9. identify unique cell cycles by ID number
            isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');              % isDrop; 1 = birth, 0 when not
            curveFinder = getGrowthParameter(conditionData_fullOnly,'curveFinder');    % curveFinder, ID of full cell cycles
            unique_cc = curveFinder(isDrop == 1);
            
            
            
            % 10. remove birth times earlier than 3 h
            timestamp_hr = getGrowthParameter(conditionData_fullOnly,'timestamp')/3600;   % raw timestamp (sec), convert to h
            birthTimes = timestamp_hr(isDrop == 1);
            
            birthTimes_final = birthTimes(birthTimes > 3);
            unique_cc_final = unique_cc(birthTimes > 3);
            clear unique_cc unique_cc_post3 birthTimes birthTimes_post3 isDrop
            
            
            
            % 10. collect tau and lambda for each cell cycle
            interdiv_fullOnly = getGrowthParameter(conditionData_fullOnly,'curveDurations')/60;  % sec, convert to min
            
            curr_tau = nan(length(unique_cc_final),1);
            curr_lambda = nan(length(unique_cc_final),1);
            curr_cclengths = nan(length(unique_cc_final),1); % to later trim cell cycles of less than 5 points
            curr_tdiv = nan(length(unique_cc_final),1); % to later trim if dividing after bubbletime
            
            for cc = 1:length(unique_cc_final)
                
                currentMus = growthRates_fullOnly(curveFinder == unique_cc_final(cc));
                currentTau = interdiv_fullOnly(curveFinder == unique_cc_final(cc));
                currentTimes = timestamp_hr(curveFinder == unique_cc_final(cc));
                
                curr_cclengths(cc,1) = length(currentMus);
                curr_lambda(cc,1) = nanmean(currentMus);
                curr_tau(cc,1) = currentTau(1);
                curr_tdiv(cc,1) = currentTimes(end);
                
            end
            
            cc_lengths = [cc_lengths; curr_cclengths];
            lambda = [lambda; curr_lambda];
            tau = [tau; curr_tau];
            t_division = [t_division; curr_tdiv];
            clear cc currentTime currentVolumes curveFinder unique_cc_final
            clear growthRates_fullOnly interdiv_fullOnly
            clear curr_cclengths curr_lambda curr_tau curr_tdiv
            
        end
        
        % 11. trim cell cycles that divide after bubbletime
        lambda_postBubbles = lambda(t_division < bubbletime(condition));
        tau_postBubbles = tau(t_division < bubbletime(condition));
        cc_lengths_postBubbles = cc_lengths(t_division < bubbletime(condition));
        
        
        % 12. trim cell cycles of less than 5 points
        lambda_6plus = lambda_postBubbles(cc_lengths_postBubbles > 5);
        tau_6plus = tau_postBubbles(cc_lengths_postBubbles > 5);
        clear tau_postBubbles lambda_postBubbles tau lambda cc_lengths cc_lengths_postBubbles t_division
        
        
        
        
        % 14. determine points outside of 3 sigma from median in tau
        tau_median = median(tau_6plus);
        tau_std = std(tau_6plus);
        
        outlier_big = find(tau_6plus > (tau_median + tau_std*3));
        outlier_small = find(tau_6plus < (tau_median - tau_std*3));
        outliers = [outlier_big; outlier_small];
        
        binary_tau = ones(length(tau_6plus),1);
        binary_tau(outliers) = 0;
        
        
        % 15. determine points outsdie of 3 sigma from media in lambda
        lambda_median = median(lambda_6plus);
        lambda_std = std(lambda_6plus);
        
        outlier_big = find(lambda_6plus > (lambda_median + lambda_std*3));
        outlier_small = find(lambda_6plus < (lambda_median - lambda_std*3));
        outliers = [outlier_big; outlier_small];
        
        binary_lambda = ones(length(lambda_6plus),1);
        binary_lambda(outliers) = 0;
        
        
        % 16. remove points outside of 3 sigma in either trait
        data_summed = binary_tau + binary_lambda;
        tau_final = tau_6plus(data_summed == 2);
        lambda_final = lambda_6plus(data_summed == 2);
        
        clear outlier_big outlier_small outliers tau_median tau_std
        clear data_summed lambda_median lambda_std binary_lambda binary_tau
        
        
        
        % 17. trim points with negative mean growth rates
        lambda_final = lambda_final(lambda_final > 0);
        tau_final = tau_final(lambda_final > 0);
        
        
        
        % 18. plot
        color = rgb(palette(condition));
        
        % division size vs. birth size
        figure(1)
        subplot(2,2,condition)
        plot(lambda_final,tau_final,'o','Color',color)
        hold on
        legend(environment(condition))
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('interdivision time (min)')
        %axis([0 max_vbirth 0 max_vdiv])
        
        
        % 15. plot with colorbar indicating density
        binSize_mu = 0.2; % 1/h
        binSize_tau = 1;  % min
        
        max_tau = 150;
        max_mu = 5;
        
        bin_tau = ceil(tau_final/binSize_tau);
        bin_lambda = ceil(lambda_final/binSize_mu);
        
        
        
        bin_mat = zeros(max_tau/binSize_tau,max_mu/binSize_mu);
        test = sub2ind(size(bin_mat),bin_tau,bin_lambda);
        
        
        bins_unique = unique(test);
        out = [bins_unique,histc(test(:),bins_unique)]; % outputs linear index in column 1, number of counts per index in column 2
        
        idx = out(:,1);
        counts = out(:,2);
        
        
        bin_mat(out(:,1)) = out(:,2);
        [i,j] = ind2sub(size(bin_mat),idx);
        
        binned_lambdas = j.*binSize_mu;
        binned_taus = i.*binSize_tau;
        
        
        
        figure(2)
        subplot(2,2,condition)
        scatter(binned_lambdas,binned_taus,60,counts,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 max_mu 0 max_tau])
        title(condition)
        xlabel('mean growth rate (1/h)')
        ylabel('interdivision time (min)')
        
        
    end
    
    % 17. save plots in active folder
    %cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    figure(1)
    plotName = strcat('B3-fig1-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    plotName = strcat('B3-fig2-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
    
    % 18. save data
    
end

