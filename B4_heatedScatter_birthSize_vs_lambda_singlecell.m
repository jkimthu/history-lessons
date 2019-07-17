%% Figure B4: birth size vs lambda, heated scatter


%  Goal: plot a heated scatter of single-cell interdivision time vs mean growth rate.
%        consider only cell cycles born after 3 hrs


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, birth lengths and mean growth rates
%       d) plot raw condition data and heatmap version, for each cycle
%       e) save both figures per experiment



%  Last edit: jen, 2019 July 17
%  Commit: first commit, build on B3 but also store data



%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
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
        birthVol = [];
        birthLength = [];
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
            experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
            cd(experimentFolder)
            filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
            load(filename,'D5','T');
            
            
            % 3. compile experiment data matrix
            xy_start = xys(condition,1);
            xy_end = xys(condition,end);
            conditionData = buildDM(D5,T,xy_start,xy_end,index,expType);
            clear D5 T xy_start xy_end e
            
            
            
            % 4. isolate volume (Va), timestamp, drop, curve, and trackNum data
            volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
            timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
            isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
            curveID = getGrowthParameter(conditionData,'curveFinder');       % curve finder (ID of curve in condition)
            trackNum = getGrowthParameter(conditionData,'trackNum');         % track number (not ID from particle tracking)
            
            
            
            % 5. calculate growth rate
            growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveID,trackNum);
            growthRates = growthRates_all(:,specificColumn);
            clear volumes timestamps_sec isDrop trackNum
            
            
            
            % 6. trim data to full cell cycles ONLY
            conditionData_fullOnly = conditionData(curveID > 0,:);
            growthRates_fullOnly = growthRates(curveID > 0,:);
            clear curveID conditionData growthRates growthRates_all
            
            
            
            % 7. identify unique cell cycles by ID number
            isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');              % isDrop; 1 = birth, 0 when not
            curveFinder = getGrowthParameter(conditionData_fullOnly,'curveFinder');    % curveFinder, ID of full cell cycles
            unique_cc = curveFinder(isDrop == 1);
            
            
            
            % 8. remove birth times earlier than 3 h
            timestamp_hr = getGrowthParameter(conditionData_fullOnly,'timestamp')/3600;   % raw timestamp (sec), convert to h
            birthTimes = timestamp_hr(isDrop == 1);
            
            birthTimes_final = birthTimes(birthTimes > 3);
            unique_cc_final = unique_cc(birthTimes > 3);
            clear unique_cc unique_cc_post3 birthTimes birthTimes_post3 isDrop
            
            
            
            % 9. collect birth size and lambda for each cell cycle
            vol_fullOnly = getGrowthParameter(conditionData_fullOnly,'volume');  
            len_fullOnly = getGrowthParameter(conditionData_fullOnly,'length');
            
            curr_birthVol = nan(length(unique_cc_final),1);
            curr_birthLen = nan(length(unique_cc_final),1);
            curr_lambda = nan(length(unique_cc_final),1);
            curr_cclengths = nan(length(unique_cc_final),1); % to later trim cell cycles of less than 5 points
            curr_tdiv = nan(length(unique_cc_final),1); % to later trim if dividing after bubbletime
            
            for cc = 1:length(unique_cc_final)
                
                currentMus = growthRates_fullOnly(curveFinder == unique_cc_final(cc));
                currentVolumes = vol_fullOnly(curveFinder == unique_cc_final(cc));
                currentLengths = len_fullOnly(curveFinder == unique_cc_final(cc));
                currentTimes = timestamp_hr(curveFinder == unique_cc_final(cc));
                
                curr_cclengths(cc,1) = length(currentMus);
                curr_lambda(cc,1) = nanmean(currentMus);
                curr_birthVol(cc,1) = currentVolumes(1);
                curr_birthLen(cc,1) = currentLengths(1);
                curr_tdiv(cc,1) = currentTimes(end);
                
            end
            
            cc_lengths = [cc_lengths; curr_cclengths];
            lambda = [lambda; curr_lambda];
            birthVol = [birthVol; curr_birthVol];
            birthLength = [birthLength; curr_birthLen];
            t_division = [t_division; curr_tdiv];
            clear cc currentTimes currentVolumes currentMus currentLengths curveFinder unique_cc_final
            clear growthRates_fullOnly vol_fullOnly len_fullOnly
            clear curr_cclengths curr_lambda curr_birthVol curr_birthLen curr_tdiv
            
        end
        
        
        % 10. trim cell cycles that divide after bubbletime
        lambda_postBubbles = lambda(t_division < bubbletime(condition));
        birthVol_postBubbles = birthVol(t_division < bubbletime(condition));
        birthLen_postBubbles = birthLength(t_division < bubbletime(condition));
        cc_lengths_postBubbles = cc_lengths(t_division < bubbletime(condition));
        
        
        
        % 11. trim cell cycles of less than 5 points
        lambda_6plus = lambda_postBubbles(cc_lengths_postBubbles > 5);
        birthVol_6plus = birthVol_postBubbles(cc_lengths_postBubbles > 5);
        birthLen_6plus = birthLen_postBubbles(cc_lengths_postBubbles > 5);
        clear birthVol_postBubbles birthLen_postBubbles lambda_postBubbles 
        clear birthVol birthLength lambda cc_lengths cc_lengths_postBubbles t_division
        
        
        
        % 12. determine points outside of 3 sigma from median in birth volume
        birthVol_median = median(birthVol_6plus);
        birthVol_std = std(birthVol_6plus);
        
        outlier_big = find(birthVol_6plus > (birthVol_median + birthVol_std*3));
        outlier_small = find(birthVol_6plus < (birthVol_median - birthVol_std*3));
        outliers = [outlier_big; outlier_small];
        
        binary_birthVol = ones(length(birthVol_6plus),1);
        binary_birthVol(outliers) = 0;
        
        

        % 13. determine points outsdie of 3 sigma from media in lambda
        lambda_median = median(lambda_6plus);
        lambda_std = std(lambda_6plus);
        
        outlier_big = find(lambda_6plus > (lambda_median + lambda_std*3));
        outlier_small = find(lambda_6plus < (lambda_median - lambda_std*3));
        outliers = [outlier_big; outlier_small];
        
        binary_lambda = ones(length(lambda_6plus),1);
        binary_lambda(outliers) = 0;
        
        
        
        % 14. remove points outside of 3 sigma in either trait
        data_summed = binary_birthVol + binary_lambda;
        birthVol_final = birthVol_6plus(data_summed == 2);
        birthLength_final = birthLen_6plus(data_summed == 2);
        lambda_final = lambda_6plus(data_summed == 2);
        
        clear outlier_big outlier_small outliers tau_median tau_std
        clear data_summed lambda_median lambda_std binary_lambda binary_tau
        
        
        
        
        % 15. trim points with negative mean growth rates
        lambda_final = lambda_final(lambda_final > 0);
        birthVol_final = birthVol_final(lambda_final > 0);
        birthLength_final = birthLength_final(lambda_final > 0);
        
        
        
        % 16. plot
        color = rgb(palette(condition));
        
        % birth volume vs. mean growth rate
        figure(1)
        subplot(2,2,condition)
        plot(lambda_final,birthVol_final,'o','Color',color)
        hold on
        legend(environment(condition))
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        %axis([0 max_vbirth 0 max_vdiv])
        
        
        % birth length vs. mean growth rate
        figure(2)
        subplot(2,2,condition)
        plot(lambda_final,birthLength_final,'o','Color',color)
        hold on
        legend(environment(condition))
        title(date)
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        
        
        % 15. plot with colorbar indicating density
        binSize_mu = 0.2; % 1/h
        binSize_vol = 0.2;  % cubic um
        binSize_length = 0.2; % um
        
        max_size = 15;
        max_mu = 5;
        
        bin_birthVol = ceil(birthVol_final/binSize_vol);
        bin_lambda = ceil(lambda_final/binSize_mu);
        
        
        bin_matrix_vol = zeros(max_size/binSize_vol,max_mu/binSize_mu);
        bin_matrix_length = zeros(max_size/binSize_length,max_mu/binSize_mu);
        
        linidx_vol = sub2ind(size(bin_matrix_vol),bin_birthVol,bin_lambda);
        linidx_len = sub2ind(size(bin_matrix_length),bin_birthVol,bin_lambda);
        
        
        bins_unique_vol = unique(linidx_vol);
        out_v = [bins_unique_vol,histc(linidx_vol(:),bins_unique_vol)]; % outputs linear index in column 1, number of counts per index in column 2
        idx_v = out_v(:,1);
        counts_v = out_v(:,2);
        
        bins_unique_len = unique(linidx_len);
        out_l = [bins_unique_len,histc(linidx_len(:),bins_unique_len)];
        idx_l = out_l(:,1);
        counts_l = out_l(:,2);
        
        
        bin_matrix_vol(out_v(:,1)) = out_v(:,2);
        [iv,jv] = ind2sub(size(bin_matrix_vol),idx_v);
        binned_lambdas_v = jv.*binSize_mu;
        binned_vol = iv.*binSize_vol;
        
        bin_matrix_length(out_l(:,1)) = out_l(:,2);
        [il,jl] = ind2sub(size(bin_matrix_length),idx_l);
        binned_lambdas_l = jl.*binSize_mu;
        binned_length = il.*binSize_length;
        
        
        
        figure(3)
        subplot(2,2,condition)
        scatter(binned_lambdas_v,binned_vol,60,counts_v,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 max_mu 0 max_size])
        title(condition)
        xlabel('mean growth rate (1/h)')
        ylabel('birth volume (cubic um)')
        
        figure(4)
        subplot(2,2,condition)
        scatter(binned_lambdas_l,binned_length,60,counts_l,'filled')
        colormap(parula) % see colormap documentation
        colorbar;
        axis([0 max_mu 0 max_size])
        title(condition)
        xlabel('mean growth rate (1/h)')
        ylabel('birth length (um)')
        
    end
    
    % 17. save plots in active folder
    %cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    figure(1)
    plotName = strcat('B4-fig1-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    plotName = strcat('B4-fig2-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
    figure(3)
    plotName = strcat('B4-fig3-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
    figure(4)
    plotName = strcat('B4-fig4-compiled-',num2str(timescale));
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
    
    % 18. store data
    birthVolumes{condition,ts} = birthVol_final;
    birthLengths{condition,ts} = birthLength_final;
    lambdas{condition,ts} = lambda_final;
    
  
end

% 19. save data
save('B4_birthSize_vs_lambda','birthVolumes','birthLengths','lambdas')