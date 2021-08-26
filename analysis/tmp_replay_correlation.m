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

    trialdur = nan(nTrls,1);
    for trl = 1:nTrls
        trialdur(trl,1) = length(data.time{trl});
    end
    
    pvals = [];
    rhos = [];
    meanamp = nan(nTrls,2,max(trialdur));
    pathfreq = [];
    pathphase = [];
    for trl = 1:nTrls

%         disp(['Trial ' num2str(trl) ' of ' num2str(nTrls) '...'])
        
        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

        [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
        
        seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
        path1 = mean(seqevidence(:,1:2),2);
        path2 = mean(seqevidence(:,3:4),2);
        x = thistrial.time{1}(1:size(seqevidence,1));
        
        if x(end) >= 5
            % ignore first second of data (items were flashing on-screen)
            xidx = x >= 1;
            path1 = path1(xidx);
            path2 = path2(xidx);
            x = x(xidx);

            [r,p] = corr(path1,path2);
            rhos(trl,1) = r;
            pvals(trl,1) = p;

            if behav.Forced(trl)==0 && ((behav.nV_1(trl)>0 && behav.nV_2(trl)<0) || (behav.nV_1(trl)<0 && behav.nV_2(trl)>0))
                if behav.nV_1(trl) > behav.nV_2(trl)
                    meanamp(trl,1,1:length(x)) = path1;
                    meanamp(trl,2,1:length(x)) = path2;
                else
                    meanamp(trl,1,1:length(x)) = path2;
                    meanamp(trl,2,1:length(x)) = path1;
                end
            end

            [freq,amp,phase] = getPhase(x,path1);
            pathfreq(trl,1) = freq(amp==max(amp));
            pathphase(trl,1) = phase(amp==max(amp));

            [freq,amp,phase] = getPhase(x,path2);
            pathfreq(trl,2) = freq(amp==max(amp));
            pathphase(trl,2) = phase(amp==max(amp));
        end
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
    thisT = behav(:,ismember(behav.Properties.VariableNames,...
        {'Practice','Block','Trial','Forced','ExpTrial','P','nV_1','nV_2','EV','Choice','Acc','RT','Subject'}));
    thisT.Replay_correlation = rhos;
    thisT.Replay_corrpvals = pvals;
    thisT.Rewarding_freq = pathfreq(:,1);
    thisT.Aversive_freq = pathfreq(:,2);
    thisT.Rewarding_phase = pathphase(:,1);
    thisT.Aversive_phase = pathphase(:,2);
    
    T = [T; thisT];
    
end

writetable(T,fullfile('D:\2020_RiskyReplay\results\timefrequency','frequency_table.csv'));

%% T-tests

vN = {'Replay_correlation','Replay_corrpvals','Rewarding_freq','Aversive_freq','Rewarding_phase','Aversive_phase','Phase_diff'};

x = [];
for s = 1:N
    
    idx = T.Subject==str2double(subjects{s}) & T.Forced==0 & T.RT >= 5;
    tmp = [table2array(T(idx,ismember(T.Properties.VariableNames,vN))) abs(rad2deg(T.Rewarding_phase(idx))-rad2deg(T.Aversive_phase(idx)))];
    x = [x; array2table(mean(tmp),'variablenames',vN)];
    
end



[H,P,CI,STATS] = ttest(x.Replay_correlation)

figure
w = .1;
[X,Y] = beeswarm(x.Replay_correlation,.05,w);
scatter(X,Y,50,'markerfacecolor','k','markeredgealpha',0,'markerfacealpha',.6); hold on
patch([0-w 0+w 0+w 0-w 0-w],[CI(2) CI(2) CI(1) CI(1) CI(2)],'w','facealpha',.3,'edgecolor','k','linewidth',1.5); hold on
plot([0-w 0+w],repmat(mean(Y),2,1),'k','linewidth',1.5); hold on
xlim([-0.5 0.5])
set(gca,'ticklength',[0 0])

cmap = [0, 237, 172 
        255, 0, 116 ]/255;

[H,P,CI,STATS] = ttest(x.Rewarding_freq,x.Aversive_freq)

figure
set(gcf,'position',[440 380 360 418])
w = .1;
for i = 1:2
    if i==1
        [X,Y] = beeswarm(x.Rewarding_freq,.15,w);
    elseif i==2
        [X,Y] = beeswarm(x.Aversive_freq,.15,w);
    end
    upper = mean(Y)+std(Y)/sqrt(length(Y));
    lower = mean(Y)-std(Y)/sqrt(length(Y));
    X = X+i;
    scatter(X,Y,50,'markerfacecolor',cmap(i,:),'markeredgealpha',0,'markerfacealpha',.6); hold on
    patch([i-w i+w i+w i-w i-w],[upper upper lower lower upper],'w','facealpha',.3,'edgecolor','k','linewidth',1.5); hold on
    plot([i-w i+w],repmat(mean(Y),2,1),'k','linewidth',1.5); hold on
end
xlim([0.5 2.5])
set(gca,'ticklength',[0 0])



[H,P,CI,STATS] = ttest(x.Phase_diff)

figure
w = .1;
[X,Y] = beeswarm(x.Phase_diff,10,w);
scatter(X,Y,50,'markerfacecolor','k','markeredgealpha',0,'markerfacealpha',.6); hold on
patch([0-w 0+w 0+w 0-w 0-w],[CI(2) CI(2) CI(1) CI(1) CI(2)],'w','facealpha',.3,'edgecolor','k','linewidth',1.5); hold on
plot([0-w 0+w],repmat(mean(Y),2,1),'k','linewidth',1.5); hold on
xlim([-0.5 0.5])
set(gca,'ticklength',[0 0])

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
