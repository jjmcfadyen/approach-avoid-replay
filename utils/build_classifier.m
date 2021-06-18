function classifier = build_classifier(filename)
% Builds a classifier using data X and labels Y saved in 'filename'
% do_build creates the classifier (true or false)
% do_test does the cross-validation (true or false)

X = [];
Y = [];

tmp = load(filename); % loads 'X' and 'Y' variables
X = tmp.X;
Y = tmp.Y;

classifier = [];

%% Settings

% Generate lambda values
gamma      = 0.05; % parameter for half-Cauchy distribution of lambdas
nLambda    = 100;  % how many lambda values to sample
lWidth     = [1e-4:1e-4:1]; % limits of sampling

rng(1);
pdf = 2 ./ (pi*gamma*(1 + (lWidth/gamma).^2)); % half-Cauchy distribution
lambdas = sort(datasample(lWidth,nLambda,'weights',pdf,'replace',false));

%% Get info from data

states = 1:6;
nStates = length(states);

nTrls = length(Y(ismember(Y,states),:));
nChan = size(X,2);

%% Build classifier

fprintf('\nBuilding classifier...')

betas = nan(nStates,nChan,nLambda);
intercepts = nan(nStates,nLambda);
parfor st = 1:nStates

    fprintf([' state ' num2str(st) ' of ' num2str(nStates) '... '])
    
    y = Y == states(st); % binarise labels

    warning('off','all')
    [B,F] = lassoglm(X, ... % X (data)
                     y, ...  % Y (labels)
                     'binomial', ...
                     'Alpha', 1,...   % relative balance between L2 and L1 penalty
                     'Lambda', lambdas, ... % how much it will downweight non-predictive features
                     'Standardize', false,...
                     'MaxIter',64);
    warning('on','all')

    betas(st,:,:) = B;
    intercepts(st,:) = F.Intercept;

end

classifier.betas = betas;
classifier.intercepts = intercepts;

end
