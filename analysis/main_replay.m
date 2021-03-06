% Main replay script

clear all
clc

%% Directories

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_meg = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';
dir_scripts = 'D:\2020_RiskyReplay\approach-avoid-replay\';
dir_classifiers = 'D:\2020_RiskyReplay\data\meg\classifiers';
dir_replay = 'D:\2020_RiskyReplay\data\meg\replay\withoutintercept';

cd(dir_scripts)

%% Settings

addpath('utils');
addpath('preprocessing')
addpath(genpath('analysis'))

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

% Downsampling
Fs = 100;

% Classification parameters
trainTimes = 0:10:300; % in ms
nT = length(trainTimes);
thisnull = 1; % proportion of data to replicate as null (zeros)
nStates = 6;

% Optimisation parameters
top_percentile = .2; % after sorting the classifiers by accuracy, within which top percentile to optimise for sequenceness?
group_best_time = 120; % from group average cross-validation accuracy

% Sequenceness parameters
U = generate_nullperms([]); 

%% Get classification accuracy & lambdas for all subjects

CV = [];
lambdas = [];
for s = 1:N

    %% Get best lambda for each classifier
    for t = 1:nT
        
        % Load cross-validation accuracy
        load(fullfile(dir_classifiers,subjects{s},...
            ['cv_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(thisnull) '.mat']));
        
        % Average over folds
        cv = squeeze(nanmean(cv));

        % Get best lambda (averaged across conditions)
        [~,best_lambda] = max(nanmean(cv));
        lambdas(s,t) = best_lambda;

        % Get minimum accuracy across all 6 states
        CV(s,t) = min(cv(:,best_lambda)); 
        
    end
end

%% Compute minimum classification accuracy for moving on to replay

q_threshold = .95; % get the top 5% of classification accuracies per participant
while true
    cv_threshold = min(quantile(CV',q_threshold));
    if any(sum(CV > cv_threshold)==N) % stop once all subjects have above-threshold accuracy for one training time
        break
    else
        q_threshold = q_threshold - .05;
    end
end

figure

subplot(1,2,1)
imagesc(CV)
caxis([(1/6) .6])
colormap([0 0 0; colours(99,'plasma')])
colorbar
title('Cross-Validated Accuracy (Minimum per State)')
xlabel('Classifier Training Time')
ylabel('Subjects')

subplot(1,2,2)
plot(trainTimes,min(CV)); hold on
plot(trainTimes,mean(CV)); hold on
plot(trainTimes,max(CV)); hold on
plot(trainTimes([1 end]),repmat(1/6,2,1),'k:'); hold on
legend({'min','mean','max','chance'})
title('Cross-Validated Accuracy (Minimum per State)')

%% Generate sequenceness for all classifier training times

for s = 1:N
    
    dir_save = fullfile(dir_replay,subjects{s});
    if ~exist(dir_save)
        mkdir(dir_save);
    end
    
    % Get task data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;
    
    if ~isfield(data,'fsample')
        data.fsample = Fs;
    end

    % Get replay for each classifier training time
    thesetimes = trainTimes(CV(s,:) > cv_threshold);
    thesetimes = setdiff(trainTimes(9:end),thesetimes);
    for t = 1:length(thesetimes)
       
%         % Create job for cluster
%         generate_jobs_replay(subjects{s},thesetimes(t),thisnull,lambdas(s,trainTimes==thesetimes(t)));
        
        disp('=========================================================')
        disp(['=== ' subjects{s} ', TRAINING TIME ' num2str(t) ' of ' num2str(length(thesetimes)) ' ======================'])
        disp('=========================================================')
        
        tic
        
        % build & save classifier
        filename = fullfile(dir_classifiers,subjects{s},['data_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']);
        classifier = build_classifier(filename);
        save(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'classifier');
        
        % get replay using best lambdas from classifier & save
        thislambda = lambdas(s,trainTimes==thesetimes(t));
        classifier.betas = squeeze(classifier.betas(:,:,thislambda));
        classifier.intercepts = classifier.intercepts(:,thislambda);
        replay = compute_replay(data,classifier,U);
        save(fullfile(dir_save,['replay_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'replay');

        toc
        
    end
end

%% Optimise classifier selection via forwards-backwards replay
% We want the largest peak, while also minimising autocorrelation

% Settings
excludeSubjects = {'396430'}; % bad classification
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

thisCV = CV(thesesubjects,:);
criteria_seqType = 'either'; % 'either' forwards or backwards, or the 'diff' fwd-bwd measure
thesetimes = trainTimes(min(thisCV) == max(min(thisCV))); % restrict the classifiers to be optimised per participant
thesetimes = [thesetimes-10 thesetimes thesetimes+10]; % look 10 ms either side

removeForced = true;
makeplots = true;

% Optimise/select
optimised_times = nan(N,1); % optimal training times
optimised_replay = [];
trialcount = nan(N,1);
for s = 1:N
   
    % check that the classifier is reliable for the selected training times
    thisT = thesetimes(CV(s,ismember(trainTimes,thesetimes)) > (1/6));
    
    if length(thesetimes) > 1
        disp(['Optimising ' num2str(thisT) ' for ' subjects{s} '...'])
    else
        disp(['Getting average replay for ' subjects{s} '...'])
    end
    
    % Get replay for each classifier training time   
    replay = [];
    for t = 1:length(thisT)
         tmp = load(fullfile(dir_replay,subjects{s},['replay_' subjects{s} '_t' num2str(thisT(t)) '_n' num2str(thisnull) '.mat']));
         replay(t,:,:,:,:) = squeeze(mean(tmp.replay,4)); % average over transitions
    end
    
    if removeForced
        
         % match up the replay data to the behavioural data
         load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
         behav = behav.task;
         
         load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
         trialinfo = merged.trialinfo;
         clear merged
         
         idx = nan(size(replay,2),1);
         for trl = 1:length(idx)
             if trialinfo(trl,1)==0
                 idx(trl,1) = find(behav.Practice==1 & behav.Trial==trialinfo(trl,2));
             else
                 idx(trl,1) = find(behav.Practice==0 & behav.Block==trialinfo(trl,1) & behav.Trial==trialinfo(trl,2));
             end
         end
         
         replay = replay(:,behav.Forced(idx)==0,:,:,:);
    end
    
    trialcount(s,1) = size(replay,2);
    
    % Plot
    if makeplots
        cmap = colours(length(thisT),'viridis');
        figure
        set(gcf,'position',[5 558 1432 420])
        for t = 1:length(thisT)
            for g = 1:3
                subplot(1,3,g)
                opts = [];
                opts.g = g;
                opts.avCol = cmap(t,:);
                opts.nullThresh = false;
                opts.subtractNull = true;
                if g==3
                    opts.subtractNull = false;
                    opts.nullThresh = false;
                end
                ss_plot(squeeze(replay(t,:,:,g,:)),opts);
            end
        end
        sgtitle([subjects{s} ' (' num2str(trialcount(s,1)) ' trials)'])
        drawnow
    end
    
    % Do optimisation if there is more than one training time
    if length(thisT) > 1
        best_traintime = optimise_replay(replay,criteria_seqType,thisT,CV(s,ismember(trainTimes,thisT)));
    else
        best_traintime = 1;
    end
    
    optimised_times(s,1) = thisT(best_traintime);
    optimised_replay(s,:,:,:) = squeeze(mean(replay(best_traintime,:,:,:,:),2)); % average over trials
    
    % Update plot
    if makeplots
        if length(thisT) > 1
            for g = 1:3
                subplot(1,3,g)
                ax = gca;
                patch([ax.XLim fliplr(ax.XLim)],sort(repmat(ax.YLim,1,2)),'w','facealpha',.75,'edgealpha',0);
                opts = [];
                opts.g = g;
                opts.avCol = cmap(best_traintime,:);
                opts.nullThresh = false;
                if g==3
                    opts.subtractNull = false;
                    opts.nullThresh = false;
                end
                ss_plot(squeeze(replay(best_traintime,:,:,g,:)),opts);
            end
            drawnow;
        end
    end
end

if length(thesetimes) > 1
    save(fullfile(dir_replay,'optimised_times.mat'),'optimised_times');
end

% Plot group average, optimised training time
excludeSubjects = {'396430'}; % bad classification
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

figure;
cc = 0;
for g = [3 1 2]
    cc = cc+1;
    subplot(1,3,cc)
    opts = [];
    opts.g = g;
%     opts.showIndividual = true;
%     opts.showOverall = false;
%     opts.indCol = colours(thisN,'viridis');
    d = squeeze(optimised_replay(thesesubjects,:,g,:));
    ss_plot(d,opts);
end

if length(thesetimes) > 1
    figure
    histogram(optimised_times);
    xlabel('Classifier Training Time (ms)')
    ylabel('Subject Count')
    title(['Optimised replay across ' num2str(thesetimes) ' ms'])
end

%% Extract the path-specific replay and put in data table

excludeSubjects = {'396430','263098','680913'}; % bad classification (1), bad performance (rest)
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

X = cell(1,N); % replay data
Y = cell(1,N); % behavioural data
for s = 1:thisN
    
    S = thesesubjects(s);
    disp(['Extracting path-specific replay for ' subjects{S} '...'])
    
    %% Get data
    % Get replay for optimised classifier training time
    load(fullfile(dir_replay,subjects{S},['replay_' subjects{S} '_t' num2str(optimised_times(s)) '_n' num2str(thisnull) '.mat']));
    
    % match up the replay data to the behavioural data
    load(fullfile(dir_behav,subjects{S},[subjects{S} '_parsedBehav.mat']))
    behav = behav.task;
    
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{S} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    trialinfo = merged.trialinfo;
    clear merged
    
    idx = nan(size(replay,1),1); % converts the [block trial] index into the row index in 'behav'
    for trl = 1:length(idx)
        if trialinfo(trl,1)==0
            idx(trl,1) = find(behav.Practice==1 & behav.Trial==trialinfo(trl,2));
        else
            idx(trl,1) = find(behav.Practice==0 & behav.Block==trialinfo(trl,1) & behav.Trial==trialinfo(trl,2));
        end
    end
    
    % crop the 'behav' variable to only include those with an associated replay trial (e.g., some MEG blocks are missing)
    behav = behav(sort(idx),:);
    
    % sort replay trials in order of behavioural trials
    [~,sortidx] = sort(idx);
    replay = replay(sortidx,:,:,:,:);
    
    % Put in Y variable
    Y{s} = behav(:,[1:3 5:6 8:10 23:26]); % just take the most relevant columns, to save space
    
    %% Organise replay
    
    % sort into good vs bad paths
    nTrls   = size(replay,1); % trials
    nPerms  = size(replay,2); % null permutations
    nSeq    = size(replay,3); % forward, backward, forward-backward
    nTrans  = size(replay,4); % number of transitions
    nLag    = size(replay,5); % no. of lags (10:10:600 ms)
    
    rewarding = nan(nTrls,nPerms,nSeq,nLag); 
    aversive = nan(nTrls,nPerms,nSeq,nLag); 
    
    thisidx = behav.nV_1 > behav.nV_2; % trials where path 1 is better than path 2
    rewarding(thisidx,:,:,:) = squeeze(mean(replay(thisidx,:,:,1:2,:),4));
    
    thisidx = behav.nV_1 < behav.nV_2; % trials where path 2 is better than path 1
    rewarding(thisidx,:,:,:) = squeeze(mean(replay(thisidx,:,:,3:4,:),4));
    
    thisidx = behav.nV_1 < behav.nV_2; % trials where path 1 is worse than path 2
    aversive(thisidx,:,:,:) = squeeze(mean(replay(thisidx,:,:,1:2,:),4));
    
    thisidx = behav.nV_1 > behav.nV_2; % trials where path 2 is worse than path 1
    aversive(thisidx,:,:,:) = squeeze(mean(replay(thisidx,:,:,3:4,:),4));
    
    % Put in X variable
    x = nan(nTrls,2,nPerms,nSeq,nLag);
    x(:,1,:,:,:) = rewarding;
    x(:,2,:,:,:) = aversive;

    X{s} = x;
    
end

% Plot rewarding vs aversive
pathDefinition = 'pos-neg'; % 'pos-neg' (only trials where one path is better than the other) 
                            % 'better-worse' (all trials, including those where both are positive or both are negative)
cmap = [0, 238, 144
        255, 0, 89]/255;
    
% % Individual subjects
% for s = 1:thisN
%     figure
%     set(gcf,'position',[680 558 1233 420])
%     cc = 0;
%     for c = 1:2
%         idx = Y{s}.Forced==0 & Y{s}.Choice==c;
%         if strcmp(pathDefinition,'pos-neg')
%             idx = idx & sum([Y{s}.nV_1 Y{s}.nV_2] > 1,2) == 1;
%         end
%         for g = 1:3
%             cc = cc+1;
%             subplot(2,3,cc)
%             for i = 1:2
%                 opts = [];
%                 opts.avCol = cmap(i,:);
%                 opts.g = g;
%                 d = squeeze(X{s}(idx,i,:,g,:));
%                 ss_plot(d,opts);
%             end
%         end
%     end
%     sgtitle([subjects{thesesubjects(s)} ': ' num2str(round(mean(Y{s}.Acc(idx))*100,2)) '%'])
% end

% Group average                            
pdata = [];
for s = 1:thisN
    
    B = Y{s};
    
    % remove blocks with poor accuracy
    outliers = zeros(size(B,1),1);

    blocks = unique([B.Practice B.Block],'rows');
    B.block_acc = nan(size(B,1),1);
    for b = 1:size(blocks,1)
        bidx = B.Practice==blocks(b,1) & B.Block==blocks(b,2);
        B.block_acc(bidx) = mean(B.Acc(bidx));
    end

    outliers(B.block_acc < .6) = 1;
    
    for c = 1:2 % choice
        idx = B.Forced==0 & B.Choice==c & ~outliers;
        if strcmp(pathDefinition,'pos-neg')
            idx = idx & sum([B.nV_1 B.nV_2] > 1,2) == 1;
        end
        pdata(s,c,:,:,:,:) = squeeze(mean(X{s}(idx,:,:,:,:))); % average over forced trials
    end
end

figure
for i = 1:2
    opts = [];
    opts.avCol = cmap(i,:);
    opts.nullThresh = false;
    cc = 0;
    for c = 1:2
        for g = [3 1 2]
            cc = cc+1;
            subplot(2,3,cc)
            opts.g = g;
            d = squeeze(pdata(:,c,i,:,g,:));
            ss_plot(d,opts);
        end
    end
end

%% Linear mixed effects modelling

arrangeLags = 'all'; % 'all', 'average', or 'max'
subtractNull = false;
g = 3;
block_acc_threshold = .55; % need to have at least this accuracy in the block

% Get replay
replay_window = 1:7; % lags to include
lags = 10:10:600; % in ms

T = [];
lmetrialcount = nan(thisN,1);
for s = 1:thisN
	
    % pick trials
    idx =  Y{s}.Forced==0;
    if strcmp(pathDefinition,'pos-neg')
        idx = idx & sum([Y{s}.nV_1 Y{s}.nV_2] > 1,2) == 1;
    end
    
    B = Y{s}(idx,:);
    
    % get differential replay (rewarding - aversive)
    if subtractNull
        d = [];
        for i = 1:2
            opts = [];
            opts.subtractNull = true;
            opts.g = 3;
            tmp = ss_plot(squeeze(X{s}(idx,i,:,g,:)),opts);
            d(:,i,:) = tmp.y;
        end
        d = tmp.y;
    else
        d = squeeze(X{s}(idx,:,1,g,replay_window)); % get just this subject at the specific lags
    end
    
    switch arrangeLags
        case 'all'
            d = squeeze(d(:,1,:)-d(:,2,:));
        case 'average'
            d = mean(d,3);
            d = squeeze(d(:,1)-d(:,2));
        case 'max'
            d = max(d,[],3);
            d = squeeze(d(:,1)-d(:,2));
    end
    
    % remove outliers
    outliers = abs(zscore(mean(d,2))) > 3;
    
    % remove blocks with poor accuracy
    blocks = unique([B.Practice B.Block],'rows');
    B.block_acc = nan(size(B,1),1);
    for b = 1:size(blocks,1)
        bidx = B.Practice==blocks(b,1) & B.Block==blocks(b,2);
        B.block_acc(bidx) = mean(B.Acc(bidx));
    end
    outliers(B.block_acc < block_acc_threshold) = 1;
    
    if any(outliers)
        disp(['(removing ' num2str(sum(outliers)) ' outliers)'])
    end
    
    d = d(~outliers,:);
    B = B(~outliers,:);
    
    lmetrialcount(s,1) = size(B,1);
    
    % combine with behavioural data
    n = size(d,1);
    thistable = [];
    for t = 1:size(d,2)
        tmp = B;
        tmp.Lag = repmat(lags(replay_window(t)),n,1);
        tmp.Replay = d(:,t);
        if t==1
            thistable = tmp;
        else
            thistable = [thistable; tmp];
        end
    end
    
    thistable.Subject = repmat(subjects(thesesubjects(s)),size(thistable,1),1);
    
    T = [T; thistable];
    
end

% Recode variables
T.Choice = T.Choice+10;
T.Choice(T.Choice==11) = 1;
T.Choice(T.Choice==12) = 0;
% T.Choice = categorical(T.Choice,[1 0],{'approach','avoid'});

T.Acc = categorical(T.Acc,[0 1],{'incorrect','correct'});
T.Subject = categorical(T.Subject,unique(T.Subject),unique(T.Subject));

% Factors that predict choice
glme = fitglme(T,'Choice~Replay+(Replay|Subject:Lag)','distribution','binomial')

