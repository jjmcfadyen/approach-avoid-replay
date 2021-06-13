function q_buildClassifiers(subject,timebin)
    
%% Directories

[sinfo, dinfo] = dir_cfg();

addpath(dinfo.tb_fieldtrip);
ft_defaults;

%% Parameters

% Set these
nullData   = 1; % 0 to 1, percentage of trials to add as zeros
trainTimes = 0:10:300; % time points to train on, in ms

gamma      = 0.05; % parameter for half-Cauchy distribution of lambdas
nLambda    = 100;  % how many lambda values to sample
lWidth     = [1e-4:1e-4:1]; % limits of sampling

% Constant
nStates    = 6;
imageNames = {'baby', 'backpack', 'bicycle', 'bowtie', 'car', 'cat',...
             'cupcake', 'hourglass', 'house', 'lamp', 'toothbrush', 'zebra'};
   
% Generate lambda values
rng(1);
pdf = 2 ./ (pi*gamma*(1 + (lWidth/gamma).^2)); % half-Cauchy distribution
lambdas = sort(datasample(lWidth,nLambda,'weights',pdf,'replace',false));
% figure
% histogram(lambdas);
            
%% Run

% Keep log of how many trials go into each image/state across subjects
imageCount = array2table(nan(1,length(imageNames)),'variablenames',imageNames); % 12 possible images
stateCount = nan(1,nStates); % 6 states

% Optimise lambda parameter for this subject
dir_save = fullfile(dinfo.data_meg_classifiers,subject);
if ~exist(dir_save)
    mkdir(dir_save);
end

% Load functional localiser data
[fl,event,~] = pp_cfg(subject,'FL');
[training,stateNames,plotdata] = getFL(subject,fl,event);

% Log trial counts
for st = 1:nStates
    imageCount(1,strfindcell(imageNames,stateNames{st})) = array2table(sum(training.L == st));
    stateCount(1,st) = sum(training.L == st);
end

% Loop over training time points and lambda parameters

fprintf('------------------------------------------------------------\n');
fprintf(['----- ' num2str(timebin) 'ms -----\n'])
fprintf('------------------------------------------------------------\n');

% subset data
sp = find(round(training.x,2) == round(timebin/1000,2));
if isempty(sp)
    error(['Could not match training time ' num2str(timebin) 'ms to training.x vector'])
end
X = squeeze(training.D(:,:,sp));

% scale
X = X ./ prctile(abs(X(:)),95);

% add null data
nL = round(nullData*size(X,1)); % get no. of points to add
X = [X; zeros(nL,size(X,2))];   % add to data
Y = [training.L; zeros(nL,1)];  % add to labels

% run across lambdas
output = q_cc_build(subject, X, Y, lambdas, true, false); % just do CV (don't save betas)

% save
fname = fullfile(dinfo.data_meg_classifiers,subject,[subject '_cc_train-' num2str(timebin) 'ms.mat']); % will save to $TMP directory
save(fname,'output');
   
end