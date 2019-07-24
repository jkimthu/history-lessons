%% Figure D: bootstrapping to determine whether lineage history influences cell cycle traits


%  Goal: calculate mean p-value and coef of var for all boostrap tests
%
%        given that each test varies with:
%                   (1) condition, and
%                   (2) the number of cell cycles per subset...
%
%        each of the following parts to this script assess each group
%        individually, building a final matrix by which to plot and store
%        the results of this analysis.
    


%  Traits of interest: 
%
%       a) interdivision time
%       b) birth size
%       c) V_div/V_birth ratio
%       d) mean growth rate of cell cycle
%       e) nutrient score




%  Last edit: jen, 2019 July 24
%  Commit: first commit, calculating average p-values from save replicates of bootstrapping test

%  OK let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/D_bootstrapping/')
trait_name = {'tau','Vbirth','ratio','mu','nScore','addedVolume'};
dates = {'2018-01-29';'2018-01-31';'2018-02-01'};

%% fluctuations, for length 5


for ee = 1:3

    date = dates{ee};
    load(strcat('D-',date,'-std-c1-length5.mat'))
    
    cv = lineage_stds./lineage_means;
    
    mean_cv(ee,:) = nanmean(cv);
    mean_pval(ee,:) = nanmean(pVals);
    min_pval(ee,:) = min(pVals);
    
    nlin = size(lineage_stds);
    numLineages(ee,:) = nlin(1);

end
min_pval
clear nlin cv ee date lineage_stds lineage_means pVals

%% fluctuations, for length 4


for ee = 1:3

    date = dates{ee};
    load(strcat('D-',date,'-std-c1-length4.mat'))
    
    cv = lineage_stds./lineage_means;
    
    mean_cv(ee,:) = nanmean(cv);
    mean_pval(ee,:) = nanmean(pVals);
    min_pval(ee,:) = min(pVals);
    
    nlin = size(lineage_stds);
    numLineages(ee,:) = nlin(1);

end
min_pval
clear nlin cv ee date lineage_stds lineage_means pVals

%% fluctuations, for length 3

for ee = 1:3

    date = dates{ee};
    load(strcat('D-',date,'-std-c1-length3.mat'))
    
    cv = lineage_stds./lineage_means;
    
    mean_cv(ee,:) = nanmean(cv);
    mean_pval(ee,:) = nanmean(pVals);
    min_pval(ee,:) = min(pVals);
    
    nlin = size(lineage_stds);
    numLineages(ee,:) = nlin(1);

end
min_pval
clear nlin cv ee date lineage_stds lineage_means pVals

%% fluctuations, for length 2

for ee = 1:3

    date = dates{ee};
    load(strcat('D-',date,'-std-c1-length2.mat'))
    
    cv = lineage_stds./lineage_means;
    
    mean_cv(ee,:) = nanmean(cv);
    mean_pval(ee,:) = nanmean(pVals);
    min_pval(ee,:) = min(pVals);
    
    nlin = size(lineage_stds);
    numLineages(ee,:) = nlin(1);

end
min_pval
clear nlin cv ee date lineage_stds lineage_means pVals


%%
save(strcat('D-',date,'-std-c1-length2'),'pVals','lineage_stds','lineage_means')







