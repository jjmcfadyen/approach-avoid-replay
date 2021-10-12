%% Behavioural modelling

clear all
clc

%% Directories & parameters

addpath('utils');

dir_data = 'D:\2020_RiskyReplay\data\behav';
dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

%% Fit models

m_plan      = []; % full, optimal planning (i.e., calculate the EV, with learning rate for value change at start of block)
m_learn     = []; % learn good vs bad path purely based on outcomes, ignore odd rule and just use the transition probabilities
m_neg       = []; % learn combinations of odd rule states based on the outcomes to figure out good vs bad path, then use transition probabilities
for s = 1:N
   
    % Read in data
    d = parse_behav(subjects{s},dir_data);
    
    % Exclude missed trials
    d = d(d.RT<30,:);
    
    % Get info
    nTrls = size(d,1);
    
    % PLAN MODEL
    pred = nan(nTrls,1);
    pred(d.Forced==0 & d.EV > 1) = 1;
    pred(d.Forced==0 & d.EV < 1) = 2;
    
    idx = ~isnan(pred);
    pred = pred(idx)-1;
    orig = d.Choice(idx)-1;
    
    beta = glmfit(pred,orig,'binomial','link','log');
    mu = glmval(beta,pred,'log');
    nll = -sum(log(poisspdf(orig,mu)));
    
    % LEARN MODEL    
    
    
end
