%% Figure A8: individual tau vs lambda, bubbled not binned


%  Goal: split cell cycles by signal experienced, plot mean with standard
%        deviation in both directions



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell interdivision time and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line



%  Last edit: jen, 2019 July 27
%  Commit: first commit of signals and phenotype consolidated on one plot


%  OK let's go!



%% Part 1. bin cell cycles in 2D by tau and mean growth rate (lambda)

% goal: heated scatter where color represents nutrient signal!

% strategy: 
%
%   0. initialize complete meta data
%   0. initialize plotting parameters
%   
%   1. loop through experiments to format data and plot
%   2. initialize experiment meta data
%   3. isolate experiment data
%   4. calculate lambda and nScores from compiled mus
%   5. remove data with growth rates > max growth rate (5 1/h)
%   6. plot unheated scatter to ensure that heated scatter plots accurately
%   7. assign each cc to bins based on tau and growth rate
%   8. sort cell cycles into 2D bins: tau vs growth rate
%   9. plot scatter, heated by cell cycle counts
%  10. plot scatter, heated by mean nScore
%  11. plot scatter, heated by signal type
%  12. save plots



clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A6_data.mat')
condition = 1;


% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','ForestGreen','Orange','GoldenRod'};
shape = 'o';


% 0 . initialize plot limits
plot_tau = 90;
plot_lambda = 3.5;
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
    
    
    % 4. calculate lambda and nScores from compiled mus
    
    
    eeLambdas = cellfun(@nanmean,eeMus);
    eeNscores = cellfun(@nanmean,eeSignals);
    clear eeMus
    
    
    % 5. remove data with too large growth rates
    taus = eeTaus(eeLambdas <= max_lambda);
    nSignals = eeSignals(eeLambdas <= max_lambda);
    nScores = eeNscores(eeLambdas <= max_lambda);
    lambdas = eeLambdas(eeLambdas <= max_lambda);
    clear eeNscores eeLambdas eeSignals eeTaus
    
    
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
    lambda_binned = accumarray(signalType,lambdas,[],@(x) {x});
    
    
    % 8. calculate mean and standard error for each signal class
    tau_mean = cellfun(@mean,tau_binned);
    tau_std = cellfun(@std,tau_binned);
    tau_count = cellfun(@length,tau_binned);
    tau_sem = tau_std./sqrt(tau_count);
    
    lambda_mean = cellfun(@mean,lambda_binned);
    lambda_std = cellfun(@std,lambda_binned);
    lambda_count = cellfun(@length,lambda_binned);
    lambda_sem = lambda_std./sqrt(lambda_count);
    
    
    % 9. plot mean and standard error of each signal class
    figure(1)
    for sc = 1:length(lambda_mean)
        color = rgb(palette{sc});
        if isempty(lambda_mean(sc)) == 1
            continue
        else
            errorbar(lambda_mean(sc),tau_mean(sc),tau_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            hold on
            errorbar(lambda_mean(sc),tau_mean(sc),lambda_sem(sc),'horizontal','Color',color)
        end
    end
    clear color sc tau_mean tau_sem tau_std tau_count lambda_mean lambda_sem lambda_std lambda_count
    clear signalType taus nSignals nScores lambda_binned lambdas tau_binned
    
end
clear onlyH onlyL H2L L2H HLH LHL HLHL LHLH


% 10. plot expected tau = 1/lambda relationship
x = linspace(0.1,plot_lambda,100);
y = 60./x; % convert from h to min

figure(1)
hold on
plot(x,y,'Color',rgb('SlateGray'))
ylabel('Interdivison time')
xlabel('Mean growth rate')
axis([0 plot_lambda 0 plot_tau])


% 11. save plot
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('A8-tau-v-lambda-bubbled');
saveas(gcf,plotName,'epsc')
close(gcf)









