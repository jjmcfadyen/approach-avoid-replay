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
group_best_time = 120; % from group average cross-validation accuracy

% Sequenceness parameters
U = generate_nullperms([]); 
 
%% Compute minimum classification accuracy for moving on to replay

[~,sortidx] = sort(max(CV(:,:,2),[],2),'descend');

figure

subplot(1,2,1)
imagesc(squeeze(CV(sortidx,:,2)))
caxis([(1/6) .7])
colormap([0 0 0; colours(99,'inferno')])
colorbar
title('Mean Cross-Validated Accuracy (per State)')
xlabel('Classifier Training Time (ms)')
ylabel('Subjects')
set(gca,'XTick',[1:5:length(trainTimes)]);
set(gca,'XTickLabels',strsplit(num2str(trainTimes([1:5:length(trainTimes)]))));

subplot(1,2,2)
plot(trainTimes,min(squeeze(CV(:,:,1)))); hold on
plot(trainTimes,mean(squeeze(CV(:,:,2)))); hold on
plot(trainTimes,max(squeeze(CV(:,:,3)))); hold on
plot(trainTimes([1 end]),repmat(1/6,2,1),'k:'); hold on
title('Cross-Validated Accuracy (per State)')

[~,bestTime] = max(mean(squeeze(CV(:,:,2))),[],2); % find the training time with the highest MINIMUM state accuracy, on average across subjects
ax = gca;
plot(repmat(trainTimes(bestTime),2,1),ax.YLim,'k--');
legend({'min','mean','max','chance','best time'})

%% Generate sequenceness for all classifier training times

planningType = 'during'; % 'during' planning (main analysis) or 'post' planning (transition to outcome)
temporalModulation = true;
temporalType = 'split'; % 'continuous' means use regressor for entire time period, 'split' means compare first half with second half

for s = 1:N
 
    % Get task data
    switch planningType
        case 'during'
            dir_save = fullfile(dir_replay,subjects{s});
            load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
        case 'post'
            dir_save = fullfile([dir_replay '_postplanning'],subjects{s});
            load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_post-task_' num2str(Fs) 'Hz.mat']));
    end
    data = merged;
    clear merged;
    
    if ~exist(dir_save)
        mkdir(dir_save);
    end
    
    if ~isfield(data,'fsample')
        data.fsample = Fs;
    end
    
%     % split post-planning into first 3 seconds (transition/outcome) or last 3 seconds (outcome)
%     if strcmp(planningType,'post')
%         first = data;
%         last = data;
%         for trl = 1:length(data.trial)
%             
%             idx = first.time{trl} <= 3;
%             first.trial{trl} = first.trial{trl}(:,idx);
%             first.time{trl} = first.time{trl}(idx);
%             
%             idx = last.time{trl} >= last.time{trl}(end)-3;
%             last.trial{trl} = last.trial{trl}(:,idx);
%             last.time{trl} = last.time{trl}(idx);
%         end
%     end

    % Get replay for each classifier training time
    thesetimes = 110:10:140;
    for t = 1:length(thesetimes)
       
        % Create job for cluster
        generate_jobs_replay(subjects{s},thesetimes(t),thisnull,lambdas(s,trainTimes==thesetimes(t)));
        
%         disp('=========================================================')
%         disp(['=== ' subjects{s} ', TRAINING TIME ' num2str(t) ' of ' num2str(length(thesetimes)) ' ======================'])
%         disp('=========================================================')
%         
%         tic
%         
% %         % build & save classifier
% %         filename = fullfile(dir_classifiers,subjects{s},['data_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']);
% %         classifier = build_classifier(filename);
% %         save(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'classifier');
%         load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']));
% 
%         % get replay using best lambdas from classifier & save
%         thislambda = lambdas(s,trainTimes==thesetimes(t));
%         classifier.betas = squeeze(classifier.betas(:,:,thislambda));
%         classifier.intercepts = classifier.intercepts(:,thislambda);
%         if ~temporalModulation
%             replay = compute_replay(data,classifier,U);
%         elseif strcmp(temporalType,'continuous')
%             replay = compute_replay_tm(data,classifier,U); 
%         elseif strcmp(temporalType,'split')
%             replay = compute_replay_split(data,classifier,U); 
%         end
%         
%         if strcmp(planningType,'post')
%             replay_all = replay;
%             if ~temporalModulation
%                 replay_first = compute_replay(first,classifier,U);
%                 replay_last = compute_replay(last,classifier,U);
%             else
%                 replay_first = compute_replay_tm(first,classifier,U); 
%                 replay_last = compute_replay_tm(last,classifier,U); 
%             end
%             replay = nan(3,size(replay_all,1),size(replay_all,2),size(replay_all,3),size(replay_all,4),size(replay_all,5));
%             replay(1,:,:,:,:,:) = replay_all;
%             replay(2,:,:,:,:,:) = replay_first;
%             replay(3,:,:,:,:,:) = replay_last;
%         end
%         
%         timetag = '';
%         if temporalModulation
%             timetag = '_time-modulated';
%             if strcmp(temporalType,'split')
%                 timetag = '_time-split';
%             end
%         end
%         
%         save(fullfile(dir_save,['replay' timetag '_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'replay');
% 
%         toc
        
    end
end

%% Optimise classifier selection via forwards-backwards replay
% We want the largest peak, while also minimising autocorrelation

% Settings
criteria_seqType = 'diff'; % 'either' forwards or backwards, or the 'diff' fwd-bwd measure
thesetimes = trainTimes(mean(squeeze(CV(:,:,2)))==max(mean(squeeze(CV(:,:,2))))); % find best classifier training time on average
thesetimes = [thesetimes-10 thesetimes thesetimes+10]; % look 10 ms either side
% [~,thesetimes] = max(min(squeeze(CV(:,:,1))));
% thesetimes = trainTimes(thesetimes);
% thesetimes = 130;

removeForced = true;
makeplots = false;

% Optimise/select
optimised_times = nan(N,1); % optimal training times
optimised_replay = [];
trialcount = nan(N,1);
for s = 1:N
   
    % check that the classifier is reliable for the selected training times
    if length(thesetimes) > 1
        thisT = thesetimes(CV(s,ismember(trainTimes,thesetimes),1) > (1/6));
    else
        thisT = thesetimes;
    end
    
    if length(thesetimes) > 1
        disp(['Optimising ' num2str(thisT) ' for ' subjects{s} '...'])
    else
        disp(['Getting replay for ' subjects{s} '...'])
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
excludeSubjects = {}; %{'396430'}; % bad classification
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

figure;
cc = 0;
for g = 1:3
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

%% Get reactivation

% reactivation = [];
% for s = 1:N
%     
%     disp(['Getting state reactivation for ' subjects{s} '...'])
%     
%     dir_save = fullfile(dir_replay,subjects{s});
%     if ~exist(dir_save)
%         mkdir(dir_save);
%     end
%     
%     % Get task data
%     load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
%     data = merged;
%     clear merged;
%     nTrls = length(data.trial);
%     
%     if ~isfield(data,'fsample')
%         data.fsample = Fs;
%     end
%     
%     % Get classifier
%     load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']));
% 
%     % match up the replay data to the behavioural data
%     load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
%     behav = behav.task;
%     
%     idx = nan(nTrls,1); % converts the [block trial] index into the row index in 'behav'
%     for trl = 1:length(idx)
%         if data.trialinfo(trl,1)==0
%             idx(trl,1) = find(behav.Practice==1 & behav.Trial==data.trialinfo(trl,2));
%         else
%             idx(trl,1) = find(behav.Practice==0 & behav.Block==data.trialinfo(trl,1) & behav.Trial==data.trialinfo(trl,2));
%         end
%     end
%     
%     % crop the 'behav' variable to only include those with an associated replay trial (e.g., some MEG blocks are missing)
%     behav = behav(sort(idx),:);
%     
%     % sort replay trials in order of behavioural trials
%     [~,sortidx] = sort(idx);
%     
%     % get best lambdas from classifier
%     thislambda = lambdas(s,trainTimes==thesetimes(t));
%     classifier.betas = squeeze(classifier.betas(:,:,thislambda));
%     classifier.intercepts = classifier.intercepts(:,thislambda);
%     
%     % Get predicted state reactivation per trial
%     thisreactivation = nan(nTrls,2); % rewarding, aversive paths
%     parfor trl = 1:nTrls
%    
%         C = classifier;
%         thistrl = sortidx(trl);
% 
%         % Scale data
%         X = data.trial{thistrl}'; % channels x samples
%         X = X ./ prctile(abs(X(:)),95);
% 
%         if any(isnan(X(:)))
%             idx = ~any(isnan(X));
%             X = X(:,idx);
%             C.betas = C.betas(:,idx);
%         end
% 
%         nSamples = size(X,1);
% 
%         % Apply classifier to get predicted data
% %         Y = normr(1 ./ (1 + exp(-(X*C.betas' + repmat(C.intercepts', [size(X,1) 1])))));
%         Y = X*C.betas';
% %         Y = normr(1 ./ (1 + exp(-(X*C.betas'))));
%         
%         % Get evidence that one path's reactivation is greater than the other's
%         dm = zeros(size(Y,1),6);
%         if behav.nV_1(trl) > behav.nV_2
%             dm(:,1:3) = 1;
%             dm(:,4:6) = -1;
%         else
%             dm(:,4:6) = 1;
%             dm(:,1:3) = -1;
%         end
%         
%         betas = pinv([dm ones(size(dm,1),1)])*Y;
%         
%         % Get overall activation
%         if behav.nV_1(trl) > behav.nV_2(trl)
%             rewarding   = mean(betas(1,1:3));
%             aversive    = mean(betas(1,4:6));
%         else
%             aversive    = mean(betas(1,1:3));
%             rewarding   = mean(betas(1,4:6));
%         end
%         thisreactivation(trl,:) = [rewarding aversive];
%         
%     end
%     
%     % save to overall
%     if size(behav,2)==30
%         tmp = behav(:,[1:3 5:6 8:10 23:26 30]);
%     elseif size(behav,2)==31
%         tmp = behav(:,[1:3 5:6 8:10 23:26 31]);
%     end
%     
%     tmp_rew = tmp;
%     tmp_rew.Reactivation_type = repmat({'rewarding'},nTrls,1);
%     tmp_rew.Reactivation = thisreactivation(:,1);
%     
%     tmp_avers = tmp;
%     tmp_avers.Reactivation_type = repmat({'aversive'},nTrls,1);
%     tmp_avers.Reactivation = thisreactivation(:,2);
%     
%     reactivation = [reactivation; tmp_rew; tmp_avers];
%     
% end
% 
% reactivation = reactivation(reactivation.Forced==0,:);
% 
% reactivation.Choice = reactivation.Choice+10;
% reactivation.Choice(reactivation.Choice==11) = 1; % approach
% reactivation.Choice(reactivation.Choice==12) = 0; % avoid
% 
% lme = fitglme(reactivation,'Reactivation ~ Reactivation_type*Choice + (1|Subject)');
% 
% % Plot
% plotdata = nan(s,2,2);
% for s = 1:N
%     for c = 1:2 % avoid, approach
%         idx = contains(reactivation.Subject,subjects{s}) & reactivation.Choice==c-1;
%         plotdata(s,c,1) = mean(reactivation.Reactivation(idx & contains(reactivation.Reactivation_type,'rewarding')));
%         plotdata(s,c,2) = mean(reactivation.Reactivation(idx & contains(reactivation.Reactivation_type,'aversive')));
%     end
% end
% 
% X = [];
% Y = [];
% cc=0;
% for c = 1:2 % avoid, approach
%     for p = 1:2 % rewarding, aversive
%         cc=cc+1;
%         [x,y] = beeswarm(squeeze(plotdata(:,c,p)),.03,.2);
%         X = [X x+cc];
%         Y = [Y y];
%     end
% end
% 
% figure
% cmap = [255, 182, 54
%         169, 116, 255]/255;
% cmap = repmat(cmap,2,1);
% 
% % (subject lines)
% for s = 1:N
%     plot(X(s,1:2),Y(s,1:2),'color',[0 0 0 .25]); hold on 
%     plot(X(s,3:4),Y(s,3:4),'color',[0 0 0 .25]); hold on 
% end
% 
% % (subject dots)
% for c = 1:size(Y,2)
%     scatter(X(:,c),Y(:,c),30,'markerfacecolor',cmap(c,:),...
%         'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
% end
% 
% % (boxplot)
% width = 0.25;
% plot([1+width/2 2-width/2],[median(Y(:,1)) median(Y(:,2))],...
%     'k','linewidth',1.6); hold on
% plot([3+width/2 4-width/2],[median(Y(:,3)) median(Y(:,4))],...
%     'k','linewidth',1.6); hold on
% for c = 1:size(Y,2)
%     m = median(Y(:,c));
%     q = quantile(Y(:,c),[.25 .75]);
%     e = quantile(Y(:,c),[.05 .95]);
%     patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
%         'w','facealpha',.75,'edgecolor','k'); hold on
%     plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
% end


%% Plot example of replay heatmap

% Pick a subject
s = 1;
T = 120;
maxLag = 60;

% Get task data
load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
data = merged;
clear merged

% Get classifier
load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(optimised_times(s,1)) '_n' num2str(thisnull) '.mat']));
L = lambdas(s,trainTimes==optimised_times(s,1));
B = classifier.betas(:,:,L)';
I = classifier.intercepts(:,L)';

% Find trial with best replay
load(fullfile(dir_replay,subjects{s},['replay_' subjects{s} '_t' num2str(T) '_n' num2str(thisnull) '.mat']));
[~,trls] = sort(squeeze(sum(replay(:,1,3,:,maxLag),4)));

% Get predicted data
trl = trls(end-5);
d = data.trial{trl}';
pred = d*B; %normr(1 ./ (1 + exp(-(d*B + repmat(I, [size(d,1) 1]))))); 
for st = 1:6
    pred(:,st) = smooth(pred(:,st),3);
end
figure; 
% set(gcf,'position',[4 378 1914 420])
imagesc(pred'); colormap(colours(100,'inferno'))
% caxis([1/6 max(caxis)])
xlim([0 0.4]*100)
set(gca,'ticklength',[0 0])

%% Extract the path-specific replay and put in data table

epochtype = 'planning'; % 'planning' or 'transition'

load(fullfile(dir_replay,'optimised_times.mat'));

excludeSubjects = {'263098','680913'}; %{'396430','263098','680913'}; % bad classification (1), bad performance (rest)
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

X = cell(1,thisN); % replay data
Y = cell(1,thisN); % behavioural data
for s = 1:thisN
    
    S = thesesubjects(s);
    disp(['Extracting path-specific replay for ' subjects{S} '...'])
    
    %% Get data
    % Get replay for optimised classifier training time
    switch epochtype
        case 'planning'
            load(fullfile(dir_replay,subjects{S},['replay_' subjects{S} '_t' num2str(optimised_times(S)) '_n' num2str(thisnull) '.mat']));
        case 'transition'
            load(fullfile([dir_replay '_postplanning'],subjects{S},['replay_' subjects{S} '_t' num2str(optimised_times(S)) '_n' num2str(thisnull) '.mat']));
            replay = squeeze(replay(2,:,:,:,:,:)); % DIM 1: 1=entire epoch, 2=transition/safe, 3=outcome/safe
    end
    
    % load behaviour
    load(fullfile(dir_behav,subjects{S},[subjects{S} '_parsedBehav.mat']))
    behav = behav.task;
    
    % add path recency as a variable
    halfwaypoint = find(behav.Block==6,1,'first'); % start of second half
    halfidx = ([1:size(behav,1)]' >= halfwaypoint)+1;
    
    behav.Path1Type = cell(size(behav,1),1);
    behav.Path2Type = cell(size(behav,1),1);
    
    pathassignment = nan(2,2);
    for h = 1:2
        for p = 1:2
            pathassignment(h,p) = median(behav.Outcome(halfidx==h & behav.Transition==p)) > 0;;
        end
    end
    
    if any(sum(pathassignment,2) ~= 1) || any(sum(pathassignment)~=1)
        error('Path assignment to reward/loss does not make sense')
    end
    
    pathassignmentlabels = cell(2,2);
    pathassignmentlabels(pathassignment==0) = {'aversive'};
    pathassignmentlabels(pathassignment==1) = {'rewarding'};
    
    behav.Path1Type(1:halfwaypoint) = pathassignmentlabels(1,1);
    behav.Path2Type(1:halfwaypoint) = pathassignmentlabels(1,2);
    behav.Path1Type((halfwaypoint+1):end) = pathassignmentlabels(2,1);
    behav.Path2Type((halfwaypoint+1):end) = pathassignmentlabels(2,2);
    
    behav.RewPathRecency_block = zeros(size(behav,1),1);
    behav.LossPathRecency_block = zeros(size(behav,1),1);
    blockidx = unique([behav.Practice behav.Block],'rows');
    for b = 1:size(blockidx,1)
        idx = find(behav.Practice==blockidx(b,1) & behav.Block==blockidx(b,2));
        for trl = 2:length(idx)
            prevtrials = flipud(behav(idx(1):(idx(trl)-1),:));
            if strcmp(behav.Path1Type{trl},'rewarding')
                lastreward = find(prevtrials.Transition==1,1,'first');
                lastloss = find(prevtrials.Transition==2,1,'first');
            else
                lastreward = find(prevtrials.Transition==2,1,'first');
                lastloss = find(prevtrials.Transition==1,1,'first');
            end
            if ~isempty(lastreward)
                behav.RewPathRecency_block(idx(trl)) = lastreward;
            end
            if ~isempty(lastloss)
                behav.LossPathRecency_block(idx(trl)) = lastloss;
            end
        end
    end
    
    behav.RewPathRecency_half = zeros(size(behav,1),1);
    behav.LossPathRecency_half = zeros(size(behav,1),1);
    for h = 1:2
        idx = find(halfidx==h);
        for trl = 2:length(idx)
            prevtrials = flipud(behav(idx(1):(idx(trl)-1),:));
            if strcmp(behav.Path1Type{trl},'rewarding')
                lastreward = find(prevtrials.Transition==1,1,'first');
                lastloss = find(prevtrials.Transition==2,1,'first');
            else
                lastreward = find(prevtrials.Transition==2,1,'first');
                lastloss = find(prevtrials.Transition==1,1,'first');
            end
            if ~isempty(lastreward)
                behav.RewPathRecency_half(idx(trl)) = lastreward;
            end
            if ~isempty(lastloss)
                behav.LossPathRecency_half(idx(trl)) = lastloss;
            end
        end
    end
   
%     % add prediction error after approach trials as a variable
%     behav.OverallBlock = behav.Block;
%     behav.OverallBlock(behav.Practice==true) = 0;
%     blocks = unique(behav.OverallBlock);
%     
%     switch epochtype
%         case 'planning' % use PREVIOUS prediction error
%             behav.RewProb = behav.P;
%             behav.RewProb(behav.nV_1<behav.nV_2) = 1-behav.P(behav.nV_1<behav.nV_2);
%             behav.RewVal = max([behav.nV_1 behav.nV_2],[],2);
%             behav.LossVal = min([behav.nV_1 behav.nV_2],[],2);
% 
%             behav.PrevOutcome = nan(size(behav,1),1); % difference between outcome and EV on previous trial
%             behav.PEbyEV = nan(size(behav,1),1); % difference between outcome and EV on previous trial
%             behav.PEbyProb = nan(size(behav,1),1); % inverse of the transition and the probability
%             behav.PEbyRew = nan(size(behav,1),1); % difference between outcome and REWARDING EV
%             behav.PEbyLoss = nan(size(behav,1),1); % difference between outcome and PUNISHING EV
%             for trl = 2:size(behav,1)
%                 % check that previous trial was from same block, and previous trial was an approach trial
%                 if behav.Practice(trl)==behav.Practice(trl-1) && behav.Block(trl)==behav.Block(trl-1)
%                     if behav.Choice(trl-1)==1
%                         behav.PEbyEV(trl) = behav.Outcome(trl-1) - behav.EV(trl-1);
%                         if behav.Transition(trl-1)==1
%                             behav.PEbyProb(trl) = behav.P(trl-1);
%                         elseif behav.Transition(trl-1)==2
%                             behav.PEbyProb(trl) = 1-behav.P(trl-1);
%                         end
%                         behav.PEbyRew(trl) = behav.Outcome(trl-1) - (behav.RewVal(trl-1)*behav.RewProb(trl-1));
%                         behav.PEbyLoss(trl) = behav.Outcome(trl-1) - (behav.LossVal(trl-1)*(1-behav.RewProb(trl-1)));
%                     else
%                         behav.PEbyEV(trl) = 0;
%                         behav.PEbyProb(trl) = 0;
%                     end
%                     behav.PrevOutcome(trl) = behav.Outcome(trl-1);
%                 end
%             end
% 
%             % add prediction error with different exponential decay rates
%             cc = 0;
%             for delta = [0.1 0.2 0.3 0.5 0.7 1 1.5 3 5 10 15 25] 
%                 cc = cc+1;
%                 behav.(['PEbyEV_e' num2str(cc)]) = nan(size(behav,1),1);
%                 for block = 1:length(blocks)
% 
%                     blockidx = behav.OverallBlock==blocks(block);
% 
%                     blockPE = behav.Outcome(blockidx) - behav.EV(blockidx);
% 
%                     weightedPE = nan(length(blockPE),1);
%                     for trl = 2:length(blockPE)
%                         thisdecay = normalise(exp(delta*[1:trl]'));
%                         weightedPE(trl) = nansum(blockPE(1:trl) .* thisdecay);
%                     end
% 
%                     thisidx = find(blockidx);
%                     behav.(['PEbyEV_e' num2str(cc)])(thisidx(2:end)) = weightedPE(1:end-1);
% 
%                 end
%             end
%         case 'transition'
%             behav.PEbyEV = nan(size(behav,1),1);
%             for block = 1:length(blocks)
%                 blockidx = behav.OverallBlock==blocks(block);
%                 behav.PEbyEV(blockidx) = behav.Outcome(blockidx) - behav.EV(blockidx);
%             end
%     end
    
    % match up the replay data to the behavioural data
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
    Y{s} = behav(:,ismember(behav.Properties.VariableNames,...
        {'Practice','Block','Trial','Forced','ExpTrial','P','nV_1','nV_2','EV','Choice','RT','Acc','Outcome',...
        'Path1Type','Path2Type','RewPathRecency_block','LossPathRecency_block','RewPathRecency_half','LossPathRecency_half',...
        'PEbyEV','PEbyProb','PrevOutcome','PEbyRew','PEbyLoss'})); % just take the most relevant columns, to save space
    Y{s} = [Y{s} behav(:,contains(behav.Properties.VariableNames,'PEbyEV_e'))];
    
    % Add learning of path value
    alphas = 0.1:0.2:0.9;
    pathlearning = zeros(length(alphas),size(Y{s},1),2);
    for a = 1:length(alphas)
        for i = 2:size(Y{s},1)

            prev_value = squeeze(pathlearning(a,i-1,:))';

            prev_outcome = [0 0];
            if behav.Transition(i) > 0
                prev_outcome(behav.Transition(i)) = behav.Outcome(i);
                pathlearning(a,i,:) = alphas(a)*(prev_outcome - prev_value);
            else
                pathlearning(a,i,:) = pathlearning(a,i-1,:); 
            end

        end
        Y{s}.(['RW_a' num2str(round(alphas(a)*10)) '_path1']) = squeeze(pathlearning(a,:,1))';
        Y{s}.(['RW_a' num2str(round(alphas(a)*10)) '_path2']) = squeeze(pathlearning(a,:,2))';
    end
    
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

%% Plot rewarding vs aversive
pathDefinition = 'pos-neg'; % 'pos-neg' (only trials where one path is better than the other) 
                            % 'better-worse' (all trials, including those where both are positive or both are negative)
cmap = [0, 238, 144
        255, 0, 89
        125 125 125]/255;

% Group average                            
pdata = [];
for s = 1:thisN
    
    B = Y{s};
    R = X{s};
    
    % remove forced-choice trials
    ridx = B.Forced==0;
    B = B(ridx,:);
    R = R(ridx,:,:,:,:);
    
    % remove blocks with poor accuracy
    outliers = zeros(size(B,1),1);

%     blocks = unique([B.Practice B.Block],'rows');
%     B.block_acc = nan(size(B,1),1);
%     for b = 1:size(blocks,1)
%         bidx = B.Practice==blocks(b,1) & B.Block==blocks(b,2);
%         B.block_acc(bidx) = mean(B.Acc(bidx));
%     end
% 
%     outliers(B.block_acc < .55) = 1;
%     alloutliers(s,:) = mean(outliers);
%     if any(outliers)
%         disp(['Removing ' num2str(round(mean(outliers)*100,2)) '% trials with poor block accuracy'])
%     end
    
    for c = 1:2 % choice
        idx = B.Choice==c & ~outliers;
        if strcmp(pathDefinition,'pos-neg')
            idx = idx & sum([B.nV_1 B.nV_2] > 1,2) == 1;
        end
        pdata(s,c,:,:,:,:) = squeeze(mean(R(idx,:,:,:,:))); % average over forced trials
    end
end

figure
for i = 1:2
    opts = [];
    opts.avCol = cmap(i,:);
    opts.nullThresh = true;
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

figure
cc = 0;
for c = 1:2
    opts = [];
    opts.avCol = cmap(3,:);
    opts.nullThresh = true;
    for g = [3 1 2]
        cc = cc+1;
        subplot(2,3,cc)
        opts.g = g;
        d = squeeze(pdata(:,c,1,:,g,:)) - squeeze(pdata(:,c,2,:,g,:));
        ss_plot(d,opts);
    end
end

%% Linear mixed effects modelling

pathDefinition = 'pos-neg';
arrangeLags = 'all'; % 'all', 'average', or 'max'
subtractNull = false;
g = 3;
block_acc_threshold = 0; % need to have at least this accuracy in the block

% Get replay
replay_window = 2:20; % lags to include
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
    
    % get replay
    if subtractNull
        d = [];
        for i = 1:2
            tmp = squeeze(X{s}(idx,i,:,g,:));
            np = squeeze(mean(tmp(:,2:end,:)));
            if g==3
                abs_np = unique(round(abs(np),4),'rows');
                npthresh = quantile(max(abs_np,[],2),.975);
                d(:,i,:) = abs(squeeze(tmp(:,1,:))) > npthresh;
                posidx = squeeze(tmp(:,1,:)) > 0;
                for j = 1:size(d,1)
                    d(j,i,~posidx(j,:)) = d(j,i,~posidx(j,:)) * -1;
                end
            else
                np = unique(round(np,4),'rows');
                npthresh = quantile(max(np,[],2),.95);
                d(:,i,:) = squeeze(tmp(:,1,:)) > npthresh;
            end
        end
    else
        d = squeeze(X{s}(idx,:,1,g,replay_window)); % get just this subject at the specific lags
    end

    d_rewarding = squeeze(d(:,1,:));
    d_aversive = squeeze(d(:,2,:));
    d_differential = d_rewarding - d_aversive;

%     % remove outliers
%     outliers = abs(zscore(mean(d,2))) > 3;
    outliers = zeros(size(d,1),1);
    
%     % remove blocks with poor accuracy
%     blocks = unique([B.Practice B.Block],'rows');
%     B.block_acc = nan(size(B,1),1);
%     for b = 1:size(blocks,1)
%         bidx = B.Practice==blocks(b,1) & B.Block==blocks(b,2);
%         B.block_acc(bidx) = mean(B.Acc(bidx));
%     end
%     outliers(B.block_acc < block_acc_threshold) = 1;
%     
%     disp(['Accuracy = ' num2str(round(mean(B.Acc)*100,2)) '%'])
%     if any(outliers)
%         disp(['(removing ' num2str(sum(outliers)) ' outliers)'])
%     end
%     
%     d = d(~outliers,:);
%     B = B(~outliers,:);
    
    lmetrialcount(s,1) = size(B,1);
    
    % combine with behavioural data
    n = size(d,1);
    thistable = [];
    for t = 1:length(replay_window)
        tmp = B;
        tmp.Lag = repmat(lags(replay_window(t)),n,1);
        tmp.Replay_rewarding = d_rewarding(:,t);
        tmp.Replay_aversive = d_aversive(:,t);
        tmp.Replay_differential = d_differential(:,t);
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
T.Choice(T.Choice==11) = 1; % approach
T.Choice(T.Choice==12) = 0; % avoid
% T.Choice = categorical(T.Choice,[1 0],{'approach','avoid'});

T.Acc = categorical(T.Acc,[0 1],{'incorrect','correct'});
T.Subject = categorical(T.Subject,unique(T.Subject),unique(T.Subject));

% Factors that predict choice
glme = fitglme(T(T.Lag==60,:),'Choice~Replay_differential+(1|Subject:Lag)','distribution','binomial')

% Save
writetable(T,['D:\2020_RiskyReplay\results\replay\replay_lme_' epochtype '.csv']);