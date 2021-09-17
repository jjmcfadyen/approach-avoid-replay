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

opts = [];
opts.removeMirror = true;
opts.pathLink = true;
opts.withinSeq = true;
pathperms = generate_nullperms(opts); 
nIterations = size(pathperms,1);

T = [];
corrdata = struct('rhos',[],'pvals',[],'lags',[],'xrhos',[],'pathfreq',[],'pathphase',[]);
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

        % get sequenceness for each path (all random iterations)
        SE = [];
        for it = 1:size(pathperms,1)
            [~, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange,pathperms(it,:));
            SE(it,:,:) = seqevidence{1};
        end
        
        thisrho = nan(nIterations,1);
        thispval = nan(nIterations,1);
        thislag = nan(nIterations,1);
        thisxrho = nan(nIterations,1);
        thispathfreq = nan(nIterations,2);
        thispathphase = nan(nIterations,2);
        for it = 1:nIterations
            
            seqevidence = squeeze(SE(it,:,:));
            seqevidence = seqevidence(sum(~isnan(seqevidence),2)>0,:);
            
            x = thistrial.time{1}(1:size(seqevidence,1));

            % ignore first second of data (items were flashing on-screen)
            xidx = x >= 1;
            thisx = x(xidx);

            if x(end)>=5

                path1 = mean(seqevidence(xidx,1:2),2);
                path2 = mean(seqevidence(xidx,3:4),2);

                [r,p] = corr(path1,path2);
                thisrho(it,1) = r;
                thispval(it,1) = p;

                [C,LAGS] = xcorr(path1,path2,'coeff');
                thislag(it,1) = abs(LAGS(C==max(C))*(1/100));
                thisxrho(it,1) = max(C);

                for path = 1:2
                    if path==1
                        [freq,amp,phase] = getPhase(thisx,path1);
                    elseif path==2
                        [freq,amp,phase] = getPhase(thisx,path2);
                    end
                    thispathfreq(it,path) = freq(amp==max(amp));
                    thispathphase(it,path) = phase(amp==max(amp));
                end
            end
        end

        rhos(:,trl) = thisrho;
        pvals(:,trl) = thispval;
        lags(:,trl) = thislag;
        xrhos(:,trl) = thisxrho;
        pathfreq(:,trl,:) = thispathfreq;
        pathphase(:,trl,:) = thispathphase;

    end
    mean(rhos(1,:),2)
    quantile(mean(rhos(2:end,:),2),.025)
  
    % Add to variables
    it = 1;
    thisT = behav(:,ismember(behav.Properties.VariableNames,...
            {'Practice','Block','Trial','Forced','ExpTrial','P','nV_1','nV_2','EV','Choice','Acc','RT','Subject'}));
    thisT.Replay_correlation = rhos(it,:)';
    thisT.Replay_corrpvals = pvals(it,:)';
    thisT.Rewarding_freq = pathfreq(it,:,1)';
    thisT.Aversive_freq = pathfreq(it,:,2)';
    thisT.Rewarding_phase = pathphase(it,:,1)';
    thisT.Aversive_phase = pathphase(it,:,2)';
    thisT.CrossCorr = xrhos(it,:)';
    thisT.CrossLag = lags(it,:)';
    T = [T; thisT];
    
    corrdata(s).rhos = rhos;
    corrdata(s).pvals = pvals;
    corrdata(s).lags = lags;
    corrdata(s).xrhos = xrhos;
    corrdata(s).pathfreq = pathfreq;
    corrdata(s).pathphase = pathphase;
    
end

save('D:\2020_RiskyReplay\results\replay_correlation.mat','T','corrdata');

%% Permutation tests

excludeSubjects = {}; %{'263098','680913'};
sidx = find(~ismember(subjects,excludeSubjects));
thisN = length(sidx);

reducedPerms = generate_nullperms([]); 
% permidx = find(ismember(reducedPerms,pathperms));
permidx = 1:size(pathperms,1);

for i = 1:4

    pd = [];
    for s = 1:thisN
        thisSubj = sidx(s);
        if i==1
            pd(s,:,:) = squeeze(nanmean(corrdata(thisSubj).rhos(permidx,:),2));
        elseif i==2
            pd(s,:,:) = squeeze(nanmean(corrdata(thisSubj).lags(permidx,:),2));
        elseif i==3
            pd(s,:,:) = squeeze(nanmean(abs(corrdata(thisSubj).pathfreq(permidx,:,1)-corrdata(thisSubj).pathfreq(permidx,:,2)),2));
        elseif i==4
            pd(s,:,:) = squeeze(nanmean(rad2deg(abs(corrdata(thisSubj).pathphase(permidx,:,1)-corrdata(thisSubj).pathphase(permidx,:,2))),2));
        end
    end

    m = squeeze(nanmean(pd));
    sem = squeeze(nanstd(pd))/sqrt(N);

    npthresh = quantile(m(2:end),[.025 .975]);
    absnpthresh = repmat(quantile(abs(m(2:end)),.95),1,2) .* [-1 1];

    figure
    histogram(m); hold on
    ax = gca;
    plot([m(1) m(1)],ax.YLim,'r','linewidth',1.6);
    plot([npthresh(1) npthresh(1)],ax.YLim,'k:','linewidth',1.6);
    plot([npthresh(2) npthresh(2)],ax.YLim,'k:','linewidth',1.6);
    if i==1
        plot([absnpthresh(1) absnpthresh(1)],ax.YLim,'k--','linewidth',1.6);
    end
    if i>1
        plot([absnpthresh(2) absnpthresh(2)],ax.YLim,'k--','linewidth',1.6);
    end
    
    if i==1
        title('Rhos')
    elseif i==2
        title('Lags')
    elseif i==3
        title('Freq')
    elseif i==4
        title('Phase')
    end
    
end

%% T-tests

excludeSubjects = {'263098','680913'};

vN = {'Replay_correlation','Replay_corrpvals','Rewarding_freq','Aversive_freq','Rewarding_phase','Aversive_phase','CrossCorr','CrossLag','Phase_diff',};

x = [];
for s = 1:N
    if ~any(strcmp(excludeSubjects,subjects{s}))
        try
            idx = T.Subject==str2double(subjects{s});
        catch
            idx = contains(T.Subject,subjects{s});
        end
        idx = idx & T.Forced==0;
        tmp = [table2array(T(idx,ismember(T.Properties.VariableNames,vN))) abs(rad2deg(T.Rewarding_phase(idx))-rad2deg(T.Aversive_phase(idx)))];
        x = [x; array2table(nanmean(tmp),'variablenames',vN)];
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

md = T;

md.Subject = categorical(T.Subject,unique(T.Subject),cellstr(num2str(unique(T.Subject))));
md.CatChoice = categorical(T.Choice,[1 2],{'approach','avoid'});
md.Certainty = T.P;
md.Certainty(T.P==.1) = .9;
md.Certainty(T.P==.3) = .7;
md.Catch = (md.nV_1>0 & md.nV_2>0) | (md.nV_1<1 & md.nV_2<0);

md = md(md.PathIteration==1 & md.Forced==0 & ~md.Catch & ~isnan(md.CrossLag),:);

lme = fitglme(md,'CrossLag~CatChoice*Certainty+(1|Subject)')


certainty = unique(md.Certainty);

y = [];
for s = 1:length(subjects)
    cc = 0;
    for c = 1:2
        for p = 1:length(certainty)
            cc = cc+1;
            y(s,cc) = nanmean(md.CrossLag(md.Subject==num2str(str2double(subjects{s})) & md.Choice==c & md.Certainty==certainty(p)));
        end
    end
end

y = array2table(y,'variablenames',{'approach_uncertain','approach_mildlycertain','approach_verycertain',...
    'avoid_uncertain','avoid_mildlycertain','avoid_verycertain'});
writetable(y,'D:\2020_RiskyReplay\results\anova_replaycorrelation.csv');

%% Correlate with theta frequency

pathperms = [1 2 3 4;  % between-paths
             1 1 2 2;  % within-path (1)
             3 3 4 4]; % within-path (2)
nIterations = size(pathperms,1); 

allrho = [];
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

    % Get path reactivation across time per trial
    reactivation = [];
    reactivation.trial = cell(nIterations,nTrls);
    reactivation.time = cell(nIterations,nTrls);
    for trl = 1:nTrls
 
        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

        [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
        
        seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
        x = thistrial.time{1}(1:size(seqevidence,1));

        % ignore first second of data (items were flashing on-screen)
        for it = 1:nIterations
            
            path1 = mean(seqevidence(:,pathperms(it,1:2)),2);
            path2 = mean(seqevidence(:,pathperms(it,3:4)),2);

            if behav.nV_1(trl) > behav.nV_2(trl)
                rewpath = path1;
                losspath = path2;
            else
                rewpath = path2;
                losspath = path1;
            end

            reactivation.trial{it,trl} = [rewpath'; losspath'];
            reactivation.time{it,trl} = x;
        end
    end
    
    % Get frequency estimates across time per trial
    load(fullfile(dir_meg,['7_merged_ds-600Hz'],[subject '_task_600Hz.mat'])); % loads 'merged' variable
    highres = merged;
    clear merged;
    
    rhos = [];
    pvals = [];
    for trl = 1:nTrls
        
        disp('============================================================================')
        disp(['=== TRIAL ' num2str(trl) ' OF ' num2str(nTrls) '===============================================']);
        disp('============================================================================')
        
        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,highres); % 600 hz data for time frequency
        
        cfg = [];
        cfg.channel    = 'MEG';
        cfg.method     = 'wavelet';
        cfg.width      = 7;
        cfg.output     = 'pow';
        cfg.foi        = 1:1:150;
        cfg.toi        = thistrial.time{1}(1):0.01:thistrial.time{1}(end);
        cfg.keeptrials = 'no';
        cfg.pad        = 'nextpow2';
        tf = ft_freqanalysis(cfg, thistrial);

        for freq = 1:length(tf.freq)
            for it = 1:nIterations
                
                % get matrices
                A = reactivation.trial{it,trl};
                B = squeeze(tf.powspctrm(:,freq,:));

                Ax = reactivation.time{it,trl};
                Bx = tf.time;

                % line up
                if length(B) > length(A)

                    fidx = findMin(Ax(1),Bx):findMin(Ax(end),Bx);
                    B = B(:,fidx); % make the two vectors start and end at the same times
                    Bx = Bx(:,fidx);

                    if length(B) ~= length(A)
                        error('need to resample/interpolate to make the same length')
                    end

                else
                    error('incorrect length')
                end

                A = abs(A(1,:) - A(2,:));
                B = mean(B);

                nanidx = isnan(B);

                [rho,pval] = corrcoef(A(~nanidx),B(~nanidx));

                rhos(it,trl,freq) = rho(2);
                pvals(it,trl,freq) = pval(2);
            end
        end
    end
    
    % Remove forced-choice & catch trials
    idx = behav.Forced==0 & behav.RT>=5 & ((behav.nV_1>0 & behav.nV_2<0) | (behav.nV_1<0 & behav.nV_2>0));
%     freqrange = find(tf.freq>=4 & tf.freq<=8);
%     
%     % Compare to null 
%     y = [];
%     null = [];
%     for freq = 1:length(freqrange)
%         y(freq) = squeeze(mean(rhos(1,idx,freqrange(freq)),2));
%         null(freq) = max(abs(squeeze(mean(rhos(2:end,idx,freqrange(freq)),2))));
%     end
%     
%     x = tf.freq(freqrange);
%     figure
%     scatter(x,y); hold on
%     scatter(x,null); hold on
%     plot(x([1 end]),repmat(quantile(null,.95),2,1));
%     plot(x([1 end]),repmat(-quantile(null,.95),2,1));

    allrho(s,:,:) = squeeze(mean(rhos(:,idx,:),2));
    
end

