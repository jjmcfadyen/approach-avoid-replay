function [output] = cc_build(subject,timebin,nullData,cType)
    
%% Parameters

% Generate lambda distribution
gamma      = 0.05; % parameter for half-Cauchy distribution of lambdas
nLambda    = 100;  % how many lambda values to sample
lWidth     = [1e-4:1e-4:1]; % limits of sampling
   
rng(1);
pdf = 2 ./ (pi*gamma*(1 + (lWidth/gamma).^2)); % half-Cauchy distribution
lambdas = sort(datasample(lWidth,nLambda,'weights',pdf,'replace',false));

%% Data

% Load functional localiser data
[fl,event,~] = pp_cfg(subject,'FL');
[training,~,~] = getFL(subject,fl,event);

if strcmp(cType,'P1vsP2')
    tmp = training.L;
    training.L(tmp < 4,:) = 1; % path 1
    training.L(tmp > 3,:) = 2; % path 2
end

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

%% Build

% set up
nStates = length(unique(Y(Y~=0)));

rng(str2num(subject));

betas = nan(nStates,size(X,2),nLambda);
intercepts = nan(nStates,nLambda);
for classifier = 1:nStates

    y = Y == classifier; % binarise labels

    fprintf(['Building classifier #' num2str(classifier) '...\n'])
    
    warning('off','all')
    [B,F] = lassoglm(X, ... % X (data)
                     y, ...  % Y (labels)
                     'binomial', ... 
                     'Alpha', 1,...   % relative balance between L2 and L1 penalty
                     'Lambda', lambdas, ... % how much it will downweight non-predictive features
                     'Standardize', false,...
                     'MaxIter',64);
    warning('on','all')

    betas(classifier,:,:) = B;
    intercepts(classifier,:) = F.Intercept;

end

output.B = betas;
output.I = intercepts;
   
end
