%% Figure J: two-factor ANOVA to test for interactions between nutrient level and timescale


%  Goal: determine effect of two independent variables on a dependent variable
%        by comparing the mean between groups
%
%        Is there an interaction between the two independent variables on
%        the dependent one?



%  In this analysis:
%
%        Independent variable #1.  nutrient level (low, high)
%        Independent variable #2.  nutrient timescale (steady, fluctuating)
%
%        Dependent variable: test 1 = interdivision time;
%                            test 2 = added volume



%  Assumptions that must be fulfilled for this statistical test:
%
%        1. dependent variable should be a continuous variable
%        2. independent variables should each consist of two or more
%           categorical, independent groups
%        3. observations (of dependent variable) should be independent
%        4. no significant outliers
%        5. dependent variable should be approximately normally distributed
%           within each categorical group
%        6. population variance of each group should be about equal



%  Strategy:
%
%        Part One. Load data from I, which gathers cell cycle data from 60 min condition
%        Part Two. Test for normality in dependent variable
%        Part Three. Test for homogeneity of variance
%        Part Four. Format data for two-way ANOVA (unbalanced) and test



%  Last edit: jen, 2019 June 19
%  Commit: first commit, two-factor ANOVA with traits, nutrient level and timescale


%  OK let's go!

%% Part One. load data

clear
clc
cd('/Users/jen/Documents/StockerLab/Data_analysis/figures_ms2/')
load('I-60min-data.mat')
clear traits_50


% to loop through each trait: interdivision time(6), added volume(3)
trait_column = [6; 3];
trait = 6;


% isolate groups of observed trait data
low_steady = traits_steady{1,2}.values(:,trait); % steady-state low
low_fluc = traits_zero(:,trait);                 % only low signal (fluc)

high_steady = traits_steady{1,4}.values(:,trait); % steady-state high
high_fluc = traits_100(:,trait);                  % only high signal (fluc)


%% Part Two. Test for normality in dependent variable (trait)


figure(1)
qqplot(low_steady)
title('qqplot: steady low')

figure(2)
qqplot(low_fluc)
title('qqplot: fluctuating with only low')

figure(3)
qqplot(high_steady)
title('qqplot: steady high')

figure(4)
qqplot(high_fluc)
title('qqplot: fluctuating with only high')

figure(11)
histogram(low_steady)
title('low, steady')
ylabel('counts')

figure(12)
histogram(low_fluc,2)
title('low, fluctuating')
ylabel('counts')
ylim([0 2])

figure(13)
histogram(high_steady)
title('high, steady')
ylabel('counts')
xlim([0 1])

figure(14)
histogram(high_fluc,10)
title('high, fluctuating')
ylabel('counts')

%% Part Three. Test for homogeneity of variance

var_steady_low = var(low_steady);
var_fluc_low = var(low_fluc);
var_steady_high = var(high_steady);
var_fluc_high = var(high_fluc);


%% Part Four. Format data for two-way ANOVA (unbalanced) and test

%  for an n-way analysis of variance, order of input:
%
%  p = anovan(1, {2 3}, 'model', 4, 'varnames',{5, 6})
%
%       1. vector of dependent variable (double)
%       2. vector of independent variable #1 (char)
%       3. vector of independent variable #2 (char)
%       4. number of ways (in this case, 2)
%       5. variable name of independent variable #1 (string)
%       6. variable name of independent variable #2 (string)


% 1. compile traits from all groups
trait_compiled = [low_steady; low_fluc; high_steady; high_fluc];
ls = length(low_steady);
lf = length(low_fluc);
hs = length(high_steady);
hf = length(high_fluc);


% 2. compile classifier of independent variable #1 (nutrient level)
classifier(1:length(trait_compiled),1) = 'a';

classifier(1:ls,1) = 'l';
classifier(ls+1:ls+lf,1) = 'l';  % where l = 'low'
classifier(ls+lf+1:end) = 'h';   %       h = 'high'


% 3. compile classifier of independent variable #1 (nutrient level)
classifier_2(1:length(trait_compiled),1) = 'a';

classifier_2(1:ls,1) = 's';
classifier_2(ls+1:ls+lf,1) = 'f';       % where s = 'steady'
classifier_2(ls+lf+1:ls+lf+hs) = 's';   %       f = 'fluctuating'
classifier_2(ls+lf+hs+1:end) = 'f';


% 4. perform test! 
p = anovan(trait_compiled,{classifier classifier_2},'model',2,'varnames',{'level','timescale'})




