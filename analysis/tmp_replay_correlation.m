clear all
clc

%% Directories

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_meg = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';
dir_scripts = 'D:\2020_RiskyReplay\approach-avoid-replay\';
dir_classifiers = 'D:\2020_RiskyReplay\data\meg\classifiers';

cd(dir_scripts)

%% Settings

addpath('utils');
addpath('preprocessing')
addpath(genpath('analysis'))

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults

% Replay onsets
lagrange = 20:10:90; % lags at which to look for onsets (in ms)
trainTimes = 0:10:300; % in ms

load('D:\2020_RiskyReplay\data\meg\replay\withoutintercept\optimised_times.mat'); % get best classifier training times per participant

%% See if replay onsets are anticorrelated across trial time

pathperms = [1 2 3 4;  % between-paths
             1 1 2 2;  % within-path (1)
             3 3 4 4]; % within-path (2)
nIterations = size(pathperms,1); 

T = [];
for s = 1:N
    
    subject = subjects{s};
    
    disp('==============================')
    disp(['=== ' subject ' ==='])
    disp('==============================')

    % Load data
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subject '_task_100Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;

    nTrls = length(data.trial);

    % Load behavioural data
    load(fullfile(dir_behav,subject,[subject '_parsedBehav.mat']))
    behav = behav.task;

    % Match with MEG data
    idx = zeros(nTrls,1);
    for trl = 1:nTrls
        if data.trialinfo(trl,1) == 0
            thispractice = 1;
            thisblock = 1;
        else
            thispractice = 0;
            thisblock = data.trialinfo(trl,1);
        end
        thistrial = data.trialinfo(trl,2);
        idx(trl,1) = find(behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial);
    end
    behav = behav(idx,:);

    % Get classifier info for this subject
    [~,lambdas] = get_bestLambdas(subject,trainTimes,1);
    load(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(optimised_times(s)) '_n1.mat'])); % loads 'classifier'
    classifier.betas = squeeze(classifier.betas(:,:,lambdas(trainTimes==optimised_times(s))));
    classifier.intercepts = squeeze(classifier.intercepts(:,lambdas(trainTimes==optimised_times(s))));

    pvals = nan(nIterations,nTrls);
    rhos = nan(nIterations,nTrls);
    lags = nan(nIterations,nTrls);
    xrhos = nan(nIterations,nTrls);
    pathfreq = nan(nIterations,nTrls,2);
    pathphase = nan(nIterations,nTrls,2);
    parfor trl = 1:nTrls
 
        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

        [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
        
        seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
        x = thistrial.time{1}(1:size(seqevidence,1));

        % ignore first second of data (items were flashing on-screen)
        xidx = x >= 1;
        thisx = x(xidx);

        thisrho = nan(nIterations,1);
        thispval = nan(nIterations,1);
        thislag = nan(nIterations,1);
        thisxrho = nan(nIterations,1);
        thispathfreq = nan(nIterations,2);
        thispathphase = nan(nIterations,2);
        if x(end)>=5
            for it = 1:size(pathperms,1)

                path1 = mean(seqevidence(xidx,pathperms(it,1:2)),2);
                path2 = mean(seqevidence(xidx,pathperms(it,3:4)),2);

                [r,p] = corr(path1,path2);
                thisrho(it,1) = r;
                thispval(it,1) = p;

                [C,LAGS] = xcorr(path1,path2,'coeff');
                thislag(it,1) = abs(LAGS(C==max(C))*(1/100));
                thisxrho(it,1) = max(C);

                [freq,amp,phase] = getPhase(thisx,path2);
                thispathfreq(it,2) = freq(amp==max(amp));
                thispathphase(it,2) = phase(amp==max(amp));
            end
        end

        rhos(:,trl) = thisrho;
        pvals(:,trl) = thispval;
        lags(:,trl) = thislag;
        xrhos(:,trl) = thisxrho;
        pathfreq(:,trl,:) = thispathfreq;
        pathphase(:,trl,:) = thispathphase;

    end
    
%     % Show correlation between paths
%     figure
%     subplot(2,1,1)
%     histogram(rhos); title('Rho')
%     subplot(2,1,2)
%     histogram(pvals); title('P values')
%     sgtitle(subjects{s})
%     
%     % Show timecourse averaged across trials (display time up to shortest trial)
%     cmap = [0, 223, 115
%         245, 0, 82]/255;
%     
%     figure;
%     x = data.time{find(trialdur==min(trialdur),1,'first')};
%     for p = 1:2
%         plot(x,squeeze(nanmean(meanamp(:,p,1:min(trialdur))))','color',cmap(p,:),'linewidth',1.1); hold on
%     end
%     title(['~' num2str(round(mean(pathfreq(:)))) ' Hz'])
%     
%     legend({[num2str(round(mean(pathfreq(:,1)))) ' Hz, ' num2str(round(mean(pathphase(:,1)),2))],...
%         [num2str(round(mean(pathfreq(:,2)))) ' Hz, ' num2str(round(mean(pathphase(:,2)),2))]})
    
    % Make table
    for it = 1:nIterations
        thisT = behav(:,ismember(behav.Properties.VariableNames,...
            {'Practice','Block','Trial','Forced','ExpTrial','P','nV_1','nV_2','EV','Choice','Acc','RT','Subject'}));
        thisT.PathIteration = repmat(it,size(thisT,1),1);
        thisT.Replay_correlation = rhos(it,:)';
        thisT.Replay_corrpvals = pvals(it,:)';
        thisT.Rewarding_freq = pathfreq(it,:,1)';
        thisT.Aversive_freq = pathfreq(it,:,2)';
        thisT.Rewarding_phase = pathphase(it,:,1)';
        thisT.Aversive_phase = pathphase(it,:,2)';
        thisT.CrossCorr = xrhos(it,:)';
        thisT.CrossLag = lags(it,:)';
        T = [T; thisT];
    end
end

writetable(T,'D:\2020_RiskyReplay\results\replay_correlation.csv');

%% T-tests

excludeSubjects = {'263098','680913'};

vN = {'Replay_correlation','Replay_corrpvals','Rewarding_freq','Aversive_freq','Rewarding_phase','Aversive_phase','CrossCorr','CrossLag','Phase_diff',};

x = [];
for s = 1:N
    if ~strcmp(subjects{s},excludeSubjects)
        try
            idx = T.Subject==str2double(subjects{s});
        catch
            idx = contains(T.Subject,subjects{s});
        end
        idx = idx & T.Forced==0 & T.RT >= 5 & T.PathIteration==1;
        tmp = [table2array(T(idx,ismember(T.Properties.VariableNames,vN))) abs(rad2deg(T.Rewarding_phase(idx))-rad2deg(T.Aversive_phase(idx)))];
        x = [x; array2table(nanmean(tmp),'variablenames',vN)];
    end
end

null = cell(1,nIterations-1);
for n = 1:nIterations-1
    for s = 1:N
        if ~strcmp(subjects{s},excludeSubjects)
            try
                idx = T.Subject==str2double(subjects{s});
            catch
                idx = contains(T.Subject,subjects{s});
            end
            idx = idx & T.Forced==0 & T.RT >= 5 & T.PathIteration==n+1;
            tmp = [table2array(T(idx,ismember(T.Properties.VariableNames,vN))) abs(rad2deg(T.Rewarding_phase(idx))-rad2deg(T.Aversive_phase(idx)))];
            null{n} = [null{n}; array2table(nanmean(tmp),'variablenames',vN)];
        end
    end
end


[H,P,CI,STATS] = ttest(x.Replay_correlation)

% figure
% w = .1;
% [X,Y] = beeswarm(x.Replay_correlation,.05,w);
% scatter(X,Y,50,'markerfacecolor','k','markeredgealpha',0,'markerfacealpha',.6); hold on
% patch([0-w 0+w 0+w 0-w 0-w],[CI(2) CI(2) CI(1) CI(1) CI(2)],'w','facealpha',.3,'edgecolor','k','linewidth',1.5); hold on
% plot([0-w 0+w],repmat(mean(Y),2,1),'k','linewidth',1.5); hold on
% xlim([-0.5 0.5])
% set(gca,'ticklength',[0 0])
% title('(anti)correlation')
% 
% plot([-.5 .5],repmat(quantile(null.Replay_correlation,.025),2,1),'r:');

y = T.Replay_correlation(~isnan(T.Replay_correlation) & T.Forced==0 & T.RT >= 5 & T.PathIteration==1 & ~ismember(T.Subject,excludeSubjects));

figure

histogram(y,20,'facecolor',[0 0 0],'facealpha',.2,'normalization','pdf','edgecolor','none'); hold on
pd = pdf(fitdist(y,'normal'),linspace(-1,1,100));
plot(linspace(-1,1,100),pd,'k');

set(gca,'ticklength',[0 0])
ax = gca;

meany = x.Replay_correlation;
meannull = [];
for n = 1:length(null)
    meannull(:,n) = null{n}.Replay_correlation;
end
meannull = quantile(abs(meannull'),.95);

plot(repmat(nanmean(meany),2,1),ax.YLim,'r','linewidth',1.4);
plot(repmat(-mean(meannull),2,1),ax.YLim,'k--','linewidth',1.4);
plot(repmat(mean(meannull),2,1),ax.YLim,'k--','linewidth',1.4);
xlim([-1 1])



y = T.CrossLag(~isnan(T.CrossLag) & T.Forced==0 & T.RT >= 5 & T.PathIteration==1 & ~ismember(T.Subject,excludeSubjects));

figure

histogram(y,20,'facecolor',[0 0 0],'facealpha',.2,'normalization','pdf','edgecolor','none'); hold on
pd = pdf(fitdist(y,'normal'),linspace(0,1.2,100));
plot(linspace(0,1.2,100),pd,'k');

set(gca,'ticklength',[0 0])
ax = gca;

meany = x.CrossLag;
meannull = [];
for n = 1:length(null)
    meannull(:,n) = null{n}.CrossLag;
end
meannull = quantile(abs(meannull'),.95);

plot(repmat(nanmean(meany),2,1),ax.YLim,'r','linewidth',1.4);
plot(repmat(mean(meannull),2,1),ax.YLim,'k--','linewidth',1.4);
xlim([0,1.1])


% Find example trial
tmp = T(T.PathIteration==1 & T.RT>=5 & T.Forced==0 & ~ismember(T.Subject,excludeSubjects),:);
tmp = tmp(tmp.Replay_correlation == min(tmp.Replay_correlation),:);

s = find(contains(subjects,tmp.Subject));
trl = T(contains(T.Subject,subjects{s}) & T.PathIteration==1,:);
trl = find(trl.ExpTrial==tmp.ExpTrial);

subject = subjects{s};
    
disp('==============================')
disp(['=== ' subject ' ==='])
disp('==============================')

% Load data
load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subject '_task_100Hz.mat'])); % loads 'merged' variable
data = merged;
clear merged;

nTrls = length(data.trial);

% Load behavioural data
load(fullfile(dir_behav,subject,[subject '_parsedBehav.mat']))
behav = behav.task;

% Match with MEG data
idx = zeros(nTrls,1);
for trl = 1:nTrls
    if data.trialinfo(trl,1) == 0
        thispractice = 1;
        thisblock = 1;
    else
        thispractice = 0;
        thisblock = data.trialinfo(trl,1);
    end
    thistrial = data.trialinfo(trl,2);
    idx(trl,1) = find(behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial);
end
behav = behav(idx,:);

% Get classifier info for this subject
[~,lambdas] = get_bestLambdas(subject,trainTimes,1);
load(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(optimised_times(s)) '_n1.mat'])); % loads 'classifier'
classifier.betas = squeeze(classifier.betas(:,:,lambdas(trainTimes==optimised_times(s))));
classifier.intercepts = squeeze(classifier.intercepts(:,lambdas(trainTimes==optimised_times(s))));

y = [];

cfg = [];
cfg.trials = trl;
thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

[onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);

seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
x = thistrial.time{1}(1:size(seqevidence,1));

if behav.nV_1(trl) > behav.nV_2
    rewpath = seqevidence(:,1:2);
    losspath = seqevidence(:,3:4);
else
    losspath = seqevidence(:,1:2);
    rewpath = seqevidence(:,3:4);
end

rewpath = mean(rewpath,2);
losspath = mean(losspath,2);

srew = smooth(rewpath,5);
sloss = smooth(losspath,5);

cmap = [0, 223, 115
    245, 0, 82]/255;

figure
% plot(x,rewpath,'color',cmap(1,:),'linewidth',1.2,'linestyle',':'); hold on
% plot(x,losspath,'color',cmap(2,:),'linewidth',1.2,'linestyle',':'); hold on
plot(x,srew,'color',cmap(1,:),'linewidth',1.4); hold on
plot(x,sloss,'color',cmap(2,:),'linewidth',1.4); hold on
set(gca,'ticklength',[0 0])
xlim([0 5])
ylim([0 .4])

%% LME

T.PhaseDiff = abs(rad2deg(T.Rewarding_phase)-rad2deg(T.Aversive_phase));

% FREQUENCY OF SEQUENCENESS
longform = [];
for s = 1:N
    
    idx = contains(T.Subject,subjects{s}) & T.Forced==0 & T.RT >= 5;
    tmp = repmat(T(idx,:),2,1);
    tmp = stack(tmp,{'Rewarding_freq','Aversive_freq'},'NewDataVariableName','Freq');
%     tmp = stack(tmp,{'Rewarding_phase','Aversive_phase'},'NewDataVariableName','Phase');
    longform = [longform; tmp];
    
end

longform.Choice = longform.Choice+10;
longform.Choice(longform.Choice==11) = 1; % approach
longform.Choice(longform.Choice==12) = 0; % avoid
longform.Subject = categorical(longform.Subject,unique(longform.Subject),unique(longform.Subject));

lme = fitlme(longform,'Freq~Freq_Indicator*Choice+RT+(1|Subject)')

% PHASE
lme = fitlme(T,'PhaseDiff~Choice+RT+(1|Subject)')

%% Plot

plotdata = cell(N,2); % 1 = x axis, 2 = y axis
choiceidx = cell(N,1);
for s = 1:N
    
    subject = subjects{s};
    
    disp('==============================')
    disp(['=== ' subject ' ==='])
    disp('==============================')

    % Load data
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subject '_task_100Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;

    nTrls = length(data.trial);

    % Load behavioural data
    load(fullfile(dir_behav,subject,[subject '_parsedBehav.mat']))
    behav = behav.task;

    % Match with MEG data
    idx = zeros(nTrls,1);
    for trl = 1:nTrls
        if data.trialinfo(trl,1) == 0
            thispractice = 1;
            thisblock = 1;
        else
            thispractice = 0;
            thisblock = data.trialinfo(trl,1);
        end
        thistrial = data.trialinfo(trl,2);
        idx(trl,1) = find(behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial);
    end
    behav = behav(idx,:);

    % Get classifier info for this subject
    [~,lambdas] = get_bestLambdas(subject,trainTimes,1);
    load(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(optimised_times(s)) '_n1.mat'])); % loads 'classifier'
    classifier.betas = squeeze(classifier.betas(:,:,lambdas(trainTimes==optimised_times(s))));
    classifier.intercepts = squeeze(classifier.intercepts(:,lambdas(trainTimes==optimised_times(s))));

    trialdur = nan(nTrls,1);
    for trl = 1:nTrls
        trialdur(trl,1) = length(data.time{trl});
    end
    plotdata{s,1} = data.time{find(trialdur==max(trialdur),1,'first')};
    
    meanamp = nan(nTrls,2,max(trialdur));
    nanidx = zeros(nTrls,1);
    for trl = 1:nTrls

        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

        [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
        
        seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
        path1 = mean(seqevidence(:,1:2),2);
        path2 = mean(seqevidence(:,3:4),2);
        x = thistrial.time{1}(1:size(seqevidence,1));
        
        % if the trial is not forced-choice & not a catch trial, we'll plot it
        if behav.Forced(trl)==0 && ((behav.nV_1(trl)>0 && behav.nV_2(trl)<0) || (behav.nV_1(trl)<0 && behav.nV_2(trl)>0))
            if behav.nV_1(trl) > behav.nV_2(trl)
                meanamp(trl,1,1:length(x)) = path1;
                meanamp(trl,2,1:length(x)) = path2;
            else
                meanamp(trl,1,1:length(x)) = path2;
                meanamp(trl,2,1:length(x)) = path1;
            end
            choiceidx{s} = [choiceidx{s}; behav.Choice(trl)];
        else
            nanidx(trl,1) = 1;
        end

    end
    plotdata{s,2} = meanamp(nanidx==0,:,:);
    
end

% Show timecourse averaged across trials (display time UP TO 5 seconds)
x = [];
y = [];
for s = 1:N
    xidx = plotdata{s,1}>=0 & plotdata{s,1}<=5;
    x(s,:) = plotdata{s,1}(xidx);
    for c = 1:2
        y(s,c,:,:) = squeeze(nanmean(plotdata{s,2}(choiceidx{s}==c,:,plotdata{s,1}>=0 & plotdata{s,1}<=5)));
    end
end
x = x(1,:);

cmap = [0, 223, 115
    245, 0, 82]/255;

figure
for c = 1:2
    for p = 1:2
        m = squeeze(mean(y(:,c,p,:)));
        sem = squeeze(std(y(:,c,p,:)));
        upper = m+sem;
        lower = m-sem;
        
        subplot(1,2,c)
        patch([x fliplr(x)],[upper' fliplr(lower')],cmap(p,:),'facealpha',.2,'edgealpha',0); hold on
        plot(x,m,'color',cmap(p,:),'linewidth',1.3); hold on
    end
end


x = data.time{find(trialdur==min(trialdur),1,'first')};
for p = 1:2
    plot(x,squeeze(nanmean(meanamp(:,p,1:min(trialdur))))','color',cmap(p,:),'linewidth',1.1); hold on
end
title(['~' num2str(round(mean(pathfreq(:)))) ' Hz'])

legend({[num2str(round(mean(pathfreq(:,1)))) ' Hz, ' num2str(round(mean(pathphase(:,1)),2))],...
    [num2str(round(mean(pathfreq(:,2)))) ' Hz, ' num2str(round(mean(pathphase(:,2)),2))]})
