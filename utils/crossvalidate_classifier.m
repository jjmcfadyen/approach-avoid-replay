function [cv, avpred] = crossvalidate_classifier(filename,directory,subject,idx,nulldata)
% Builds a classifier using data X and labels Y saved in 'filename'
% do_build creates the classifier (true or false)
% do_test does the cross-validation (true or false)

X = [];
Y = [];
X1 = [];
X2 = [];
Y1 = [];
Y2 = [];

if ~isempty(filename)
    load(filename); % loads 'X' and 'Y' variables (or X1 Y1 X2 Y2 if using different train/test sets, for temporal generalisation)
else
    
    disp(['directory = ' directory])
    disp(['subject = ' subject])
    disp(['$SGE_TASK_ID = ' num2str(idx)])
    disp(['nulldata = ' num2str(nulldata)])
    filename = fullfile(directory,subject,['data_' subject '_idx' num2str(idx) '_n' num2str(nulldata) '.mat']);
    load(filename);

    tmp = strsplit(filename,'_');
    filename = [tmp{1} '_' tmp{2} '_' tmp{3} '_t' num2str(timeidx(1) + timeidx(2)/1000) '_' tmp{5}];
    disp(['FILENAME = ' filename])

end

TG = false;
if exist('X1','var')
    TG = true;
    Y = Y1; % Y1 and Y2 should be the same anyway
    X = X1; % to get channel info, etc, which should also be the same
end

cv = [];
avpred = [];

%% Settings

% Generate lambda values
gamma      = 0.05; % parameter for half-Cauchy distribution of lambdas
nLambda    = 100;  % how many lambda values to sample
lWidth     = [1e-4:1e-4:1]; % limits of sampling

rng(1);
pdf = 2 ./ (pi*gamma*(1 + (lWidth/gamma).^2)); % half-Cauchy distribution
lambdas = sort(datasample(lWidth,nLambda,'weights',pdf,'replace',false));

% Cross-validation
nFolds = 5; % integer, or 'Inf' to do as many as possible (note that any leftover trials get put in an additional fold, so specify one less than desired)
onlyOne = nFolds==Inf;

%% Get info from data

states = 1:6;
nStates = length(states);

nTrls = length(Y(ismember(Y,states),:));
nChan = size(X,2);

tmp = strsplit(filename,filesep);
thisdir = strjoin(tmp(1:end-1),filesep);

thisfile = strsplit(tmp{end},'.mat');
thisfile = strsplit(thisfile{1},'_');
thisfile = strjoin(thisfile(2:end),'_');

%% Cross-validate
    
fprintf('Determining folds...')

if nFolds==Inf
    nTrialsPerFold = 1;
else
    nTrialsPerFold = floor(nTrls/nFolds/nStates); % how many trials (of each CONDITION) are included per fold
end

% keep track of many times a trial has been included in a fold
if nTrialsPerFold ~= nTrls
    selected = array2table([ [1:nTrls]' Y(ismember(Y,states),:) zeros(nTrls,1) ],'variablenames',{'Trial','State','Count'});
    foldlog = [];
    cc = 0;
    while true
    
        thisfold = nan(nStates,nTrialsPerFold);
        for st = 1:nStates
            tmp = selected(selected.State==states(st),:);
            tmp = tmp(tmp.Count==min(tmp.Count),:);
            if size(tmp,1)<=nTrialsPerFold
                tmp = tmp.Trial';
            else
                tmp = randsample(tmp.Trial,nTrialsPerFold)';
            end
            thisfold(st,1:length(tmp)) = tmp;
            selected.Count(ismember(selected.Trial,thisfold(st,:))) = selected.Count(ismember(selected.Trial,thisfold(st,:)))+1;
        end
    
        cc = cc + 1;
    
        foldlog{cc} = thisfold;
    
        if all(selected.Count>0) || cc >= nFolds
            break
        end
    
    end
else
    foldlog = [];
    foldlog{1} = find(ismember(Y,states));
end
fprintf([' ' num2str(length(foldlog)) ' folds \n'])

% Cycle through folds
nFolds = length(foldlog);
cv = nan(nFolds,nStates,nLambda);
avpred = nan(nFolds,nStates,nStates,nLambda);
for fold = 1:nFolds

    disp(['------ fold ' num2str(fold) ' of ' num2str(nFolds) '...'])

    % Sort out left in/out data
    loIdx = foldlog{fold}(:);
    loIdx = loIdx(~isnan(loIdx));

    liIdx = setdiff(1:length(Y),loIdx);

    liY = Y(liIdx,:);
    loY = Y(loIdx,:);

    if ~TG
        liX = X(liIdx,:);
        loX = X(loIdx,:);
    else
        liX = X1(liIdx,:); % train on first data set (t1)
        loX = X2(loIdx,:); % test on second data set (t2)
    end

    thesestates = unique(loY); % some folds don't include certain states

    acc = zeros(length(thesestates),nLambda);
    PRED = zeros(length(thesestates),length(thesestates),nLambda);
    for st = 1:length(thesestates)

        % Create classifiers using left-in data
        warning('off','all')
        [B,F] = lassoglm(liX, ... % X (data)
                         liY==thesestates(st), ...  % Y (labels)
                         'binomial', ...
                         'Alpha', 1,...   % relative balance between L2 and L1 penalty
                         'Lambda', lambdas, ... % how much it will downweight non-predictive features
                         'Standardize', false,...
                         'MaxIter',64);
        warning('on','all')

        % apply classifier to left-out data
        pred = 1 ./ (1 + exp(-(loX*B + repmat(F.Intercept, [size(loX,1) 1]))));

        % sort the predictions from highest to lowest
        [~,sortidx] = sort(pred,'descend');
        sorty = loY(sortidx);

        if onlyOne
            
            % see whether the max accuracy is for this state
            thisacc = sorty(1,:)==thesestates(st);
            
            % invalidate any instances where there is more than one max prediction
            tmp = nan(1,size(pred,2));
            for i = 1:size(tmp,2)
                tmp(i) = sum(pred(:,i) == max(pred(:,i))) > 1;
            end
            
            thisacc(tmp==1) = 0;
            
            % log prediction
            thispred = pred;
            
        else
            
            % see what proportion of the highest predictions match the correct label
            thissum = sum(loY==thesestates(st));
            thisacc = mean(sorty(1:thissum,:)==thesestates(st));

            % invalidate any instances where there is more than one max prediction
            thismax = max(pred);
            tmp = nan(size(pred,1),size(pred,2));
            for i = 1:size(tmp,1)
                tmp(i,:) = pred(i,:) == thismax;
            end

            thisacc(sum(tmp) > thissum) = 0;

            % log prediction
            thispred = nan(nStates,nLambda);
            for st2 = 1:nStates
                thispred(st2,:) = mean(pred(loY==st2,:));
            end
            
        end

        % store in variable
        acc(st,:) = thisacc;
        PRED(st,:,:) = thispred;

    end

    tmpacc = squeeze(cv(fold,:,:));
    tmpacc(ismember(states,thesestates),:) = acc;
    cv(fold,:,:) = tmpacc;
    avpred(fold,:,:,:) = PRED;

end

% Display accuracy averaged over folds & states
disp(['Acc = ' num2str(round(max(squeeze(nanmean(nanmean(cv,2))))*100,2)) '%'])

%% Save output

% save(fullfile(thisdir,['cv_' thisfile '.mat']),'cv','avpred');
save(fullfile(thisdir,['cv_' thisfile '.mat']),'cv');

end
