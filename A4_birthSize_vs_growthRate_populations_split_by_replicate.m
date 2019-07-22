%% Figure A4: mean birth size vs mean growth rate


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




%  Last edit: jen, 2019 July 22
%  Commit: first commit, like A3 but plots each replicate individually to
%  get sense of deviation from best fit line


%  OK let's go!



%% Part 1. initialize analysis for each fluc replicate

% goal: plot of newborn cell size vs growth rate,
%       overlaying population and single cell data


% strategy: 
%
%       i. accumulate data from each steady condition
%          accumulate data from each replicate of each fluctuating condition
%          conditions: 3 steady conditions and 13 fluctuating replicates (16 total)
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
load('A3_data_1sigma.mat')
clear size mu bLength

% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};

shape = 'o';
sigmas = 1;

%% Part 2. assemble data per replicate

% 1. accumulate data from each condition
fluc = 1; % row number in data structure

counter = 0;
for ee = 2:5 % 2-5 correspond to fluctuating environments
    
    condition = environment_order{ee};
    
    % fluctuating environment! concatenate based on timescale
    if condition == 30
        
        idx = [2,3,4]; % ID of experiments with this fluc timescale
        
        % loop through experiments and store timescale data
        for arrayIndex = 1:length(idx)
            
            counter = counter + 1;
            expt = find(exptArray == idx(arrayIndex));
            
            % isolate data
            expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
            expt_volume = compiled_birthSize{expt,1}{fluc,1};
            expt_length = compiled_birthLength{expt,1}{1,fluc};
            
            % if data exists, calculate mean of each cell cycle
            if ~isempty(expt_mu)
                expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
            end
            
            
            % isolate cycles within 3 st dev of mean
            expt_mean = nanmean(expt_mu);
            expt_std = nanstd(expt_mu);
            
            lower = expt_mu < expt_mean + (expt_std * sigmas);
            upper = expt_mu > expt_mean - (expt_std * sigmas);
            combined = lower + upper;
            
            range_mu = expt_mu(combined == 2);
            range_volume = expt_volume(combined == 2);
            range_length = expt_length(combined == 2);
            clear expt_mean expt_std lower upper combined
            
            
            % store individual cell cycle values
            lambdas{counter} = range_mu;
            birthVolumes{counter} = range_volume;
            birthLength{counter} = range_length;
            clear expt_mu expt_volume expt_length
            
        end
        
          
    elseif condition == 300
        
        idx = [5,6,7]; % ID of experimennts with this fluc timescale
        
        % loop through experiments and store timescale data
        for arrayIndex = 1:length(idx)
            
            counter = counter + 1;
            expt = find(exptArray == idx(arrayIndex));
            
            % isolate data
            expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
            expt_volume = compiled_birthSize{expt,1}{fluc,1};
            expt_length = compiled_birthLength{expt,1}{1,fluc};
            
            % if data exists, calculate mean of each cell cycle
            if ~isempty(expt_mu)
                expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
            end
            
            
            % isolate cycles within 3 st dev of mean
            expt_mean = nanmean(expt_mu);
            expt_std = nanstd(expt_mu);
            
            lower = expt_mu < expt_mean + (expt_std * sigmas);
            upper = expt_mu > expt_mean - (expt_std * sigmas);
            combined = lower + upper;
            
            range_mu = expt_mu(combined == 2);
            range_volume = expt_volume(combined == 2);
            range_length = expt_length(combined == 2);
            clear expt_mean expt_std lower upper combined
            
            
            % store individual cell cycle values
            lambdas{counter} = range_mu;
            birthVolumes{counter} = range_volume;
            birthLength{counter} = range_length;
            clear expt_mu expt_volume expt_length
            

        end
  
        
    elseif condition == 900
        
        idx = [9,10,11,12]; % ID of experimennts with this fluc timescale

        
        % loop through experiments and store timescale data
        for arrayIndex = 1:length(idx)
            
            counter = counter + 1;
            expt = find(exptArray == idx(arrayIndex));
            
            % isolate data
            expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
            expt_volume = compiled_birthSize{expt,1}{fluc,1};
            expt_length = compiled_birthLength{expt,1}{1,fluc};
            
            % if data exists, calculate mean of each cell cycle
            if ~isempty(expt_mu)
                expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
            end
            
            % isolate cycles within 3 st dev of mean
            expt_mean = nanmean(expt_mu);
            expt_std = nanstd(expt_mu);
            
            lower = expt_mu < expt_mean + (expt_std * sigmas);
            upper = expt_mu > expt_mean - (expt_std * sigmas);
            combined = lower + upper;
            
            range_mu = expt_mu(combined == 2);
            range_volume = expt_volume(combined == 2);
            range_length = expt_length(combined == 2);
            clear expt_mean expt_std lower upper combined
            
            % store individual cell cycle values
            lambdas{counter} = range_mu;
            birthVolumes{counter} = range_volume;
            birthLength{counter} = range_length;
            clear expt_mu expt_volume expt_length
            
        end
        
        
    elseif condition == 3600
        
        idx = [13,14,15]; % ID of experimennts with this fluc timescale
        

        % loop through experiments and store timescale data
        for arrayIndex = 1:length(idx)
            
            counter = counter + 1;
            expt = find(exptArray == idx(arrayIndex));
            
            % isolate data
            expt_mu = compiled_mu{expt,1}{fluc,1}; % note: mu is all instananeous vals in each cell cycle
            expt_volume = compiled_birthSize{expt,1}{fluc,1};
            expt_length = compiled_birthLength{expt,1}{1,fluc};
            
            % if data exists, calculate mean of each cell cycle
            if ~isempty(expt_mu)
                expt_mu = cellfun(@nanmean,compiled_mu{expt,1}{fluc,1});
            end
            
            
            % isolate cycles within 3 st dev of mean
            expt_mean = nanmean(expt_mu);
            expt_std = nanstd(expt_mu);
            
            lower = expt_mu < expt_mean + (expt_std * sigmas);
            upper = expt_mu > expt_mean - (expt_std * sigmas);
            combined = lower + upper;
            
            range_mu = expt_mu(combined == 2);
            range_volume = expt_volume(combined == 2);
            range_length = expt_length(combined == 2);
            clear expt_mean expt_std lower upper combined
            
            
            % store individual cell cycle values
            lambdas{counter} = range_mu;
            birthVolumes{counter} = range_volume;
            birthLength{counter} = range_length;
            clear expt_mu expt_volume expt_length
            

        end
        
    end
end

clear fluc idx condition ee expt arrayIndex
clear counter range_mu range_volume range_length

%% Part 3. plot population means per fluc replicate with compiled steady data


% 1. calculate mean birth size and growth rate for each condition (population data)

% assemble data based on environment order
rep_birthVols = cellfun(@mean,birthVolumes);   % population birth volume of each fluctuating replicate
rep_lambdas = cellfun(@nanmean,lambdas);       % population growth rate of each fluctuating replicate
rep_birthLengths = cellfun(@mean,birthLength); % population birth length of each fluctuating replicate
% save hard earned data!
%save('A4_data_1sigma.mat','compiled_birthSize','compiled_birthLength','compiled_mu','exptArray','population_birthSize','population_mu','population_birthLength','birthVolumes','lambdas','birthLength')



% 2. plot each compiled population point
figure(1)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthSize(cc)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
    hold on
    
end
legend('low','30','300','900','3600','ave','high')
title('population growth law')
xlabel('mean lambda')
ylabel('mean birth volume')
ylim([0.5 2])
xlim([-0.1 4])
clear color cc

figure(2)
for cc = 1:length(population_mu)
    
    color = rgb(palette(cc));
    plot(population_mu(cc),log(population_birthLength(cc)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
    hold on
    
end
legend('low','30','300','900','3600','ave','high')
title('population growth law')
clear color cc



% 3. fit best line to compiled steady data

% 	i. determine which data points are steady = 1 or fluctuating = 0
for eo = 1:length(environment_order)
    isSteady(eo) = ischar(environment_order{eo});
end
clear eo

%   ii. isolate steady data
steady_mu = population_mu(isSteady == 1);
steady_birthLength = population_birthLength(isSteady == 1);
steady_birthVolume = population_birthSize(isSteady == 1);
clear isSteady

% 	iii. best fit line
fit_length = polyfit(steady_mu,log(steady_birthLength),1);
fit_vol = polyfit(steady_mu,log(steady_birthVolume),1);
clear steady_birthLength steady_birthVolume

%   iv.  plot best fit line
x = linspace(steady_mu(1),steady_mu(end),10);
yv = fit_vol(1)*x + fit_vol(2);
yl = fit_length(1)*x + fit_length(2);
clear steady_mu

figure(1)
hold on
plot(x,yv,'Color',rgb('SlateGray'))

xlabel('mean lambda')
ylabel('ln(mean birth volume)')
legend('low','30','300','900','3600','ave','high')
text(2, 1.65, strcat('y=',num2str(fit_vol(1)),'x+',num2str(fit_vol(2))))

figure(2)
hold on
plot(x,yl,'Color',rgb('SlateGray'))

xlabel('mean lambda')
ylabel('ln(mean birth length)')
legend('low','30','300','900','3600','ave','high')
text(2, 1.35, strcat('y=',num2str(fit_length(1)),'x+',num2str(fit_length(2))))

clear fit_length fit_vol x yl yv


% 4. plot each replicate of fluctuating data in filled circles
tscale = [30;30;30;300;300;300;900;900;900;900;3600;3600;3600];
for ts = 1:length(tscale)
    
    % determine color based on timescale
    if tscale(ts) == 30
        color = rgb('DarkTurquoise');
    elseif tscale(ts) == 300
        color = rgb('SteelBlue');
    elseif tscale(ts) == 900
        color = rgb('DeepSkyBlue');
    elseif tscale(ts) == 3600
        color = rgb('DodgerBlue');
    end
    
    figure(1)
    hold on
    plot(rep_lambdas(ts),log(rep_birthVols(ts)),'Color',color,'Marker',shape,'MarkerFaceColor',color,'MarkerSize',10)
    
    figure(2)
    hold on
    plot(rep_lambdas(ts),log(rep_birthLengths(ts)),'Color',color,'Marker',shape,'MarkerFaceColor',color,'MarkerSize',10)
    
end
figure(1) % vol
axis([0.7 3.3 0.6 1.8])

figure(2) % length
axis([0.7 3.3 0.7 1.5])

