function thisT = compute_timefrequency(subject,optimised_time,directories,makePlots)
% directories = structure with dir_meg, dir_behav, dir_classifiers

%% Settings

lagrange = 20:10:90; % lags at which to look for onsets (in ms)
trainTimes = 0:10:300; % in ms

%% Compute time frequency for this subject

if ~exist(directories.dir_save)
    mkdir(directories.dir_save)
end

disp('==============================')
disp(['=== ' subject ' ==='])
disp('==============================')

% Load data
load(fullfile(directories.dir_meg,['7_merged_ds-600Hz'],[subject '_task_600Hz.mat'])); % loads 'merged' variable
data = merged;
clear merged;

load(fullfile(directories.dir_meg,['7_merged_ds-100Hz'],[subject '_task_100Hz.mat'])); % loads 'merged' variable
lowdata = merged;
clear merged;

nTrls = length(data.trial);

% Load behavioural data
load(fullfile(directories.dir_behav,subject,[subject '_parsedBehav.mat']))
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
load(fullfile(directories.dir_classifiers,subject,['classifier_' subject '_t' num2str(optimised_time) '_n1.mat'])); % loads 'classifier'
classifier.betas = squeeze(classifier.betas(:,:,lambdas(trainTimes==optimised_time)));
classifier.intercepts = squeeze(classifier.intercepts(:,lambdas(trainTimes==optimised_time)));

% Do time-frequency for each trial
TF = cell(1,nTrls);
bTF = cell(1,nTrls);
rTF = cell(1,nTrls);
rbTF = cell(1,nTrls);
freq = cell(1,nTrls);
timex = cell(1,nTrls);
for trl = 1:nTrls % DO NOT DO THIS IN PARALLEL - the time-frequency computation uses LOADS of CPU (50% usage)
    
    disp(['||| TRIAL ' num2str(trl) ' of ' num2str(nTrls) ' |||'])
    
    cfg = [];
    cfg.trials = trl;
    thistrial = ft_selectdata(cfg,data);
    thistrial.trialinfo = data.trialinfo(trl,:);
    
    cfg = [];
    cfg.channel    = 'MEG';
    cfg.method     = 'wavelet';
    cfg.width      = 5;
    cfg.output     = 'pow';
    cfg.foi        = 2:1:150;
    cfg.toi        = thistrial.time{1}(1):0.01:thistrial.time{1}(end);
    cfg.keeptrials = 'no';
    cfg.pad        = 'nextpow2';
    
    tf = ft_freqanalysis(cfg, thistrial);
    TF{trl} = squeeze(nanmean(tf.powspctrm)); % average across channels
    freq{trl} = tf.freq;
    timex{trl} = tf.time;
    
    cfg = [];
    cfg.baseline = [-.3 -.1];
    btf = ft_freqbaseline(cfg,tf);
    bTF{trl} = squeeze(nanmean(btf.powspctrm)); % average across channels
    
    %{
        figure
        imagesc(TF{trl},'XData',tf.time,'YData',btf.freq);
        colormap(colours(100,'inferno'))
        view(180,90);
        set(gca,'xdir','reverse')
        hold on
        plot([0 0],tf.freq([1 end]),'w:');
        plot(repmat(thistrial.time{1}(end)-.5,2,1),tf.freq([1 end]),'w:');
    %}
 
    cfg = [];
    cfg.trials = trl;
    thislow = ft_selectdata(cfg,lowdata); % 100 hz data for replay
    thislow.trialinfo = thislow.trialinfo(trl,:);
    
    bestlag = 20:10:90;
    [onsets, seqevidence] = get_replayOnsets(thislow,classifier,bestlag);
    onsets = onsets(onsets.Path==0,:);
    otimes = zeros(1,length(btf.time));
    for o = 1:size(onsets,1)
        otimes(findMin(onsets.Onset_time(o),btf.time)) = 1;
    end
    
    % see if evidence of each path's sequenceness is anticorrelated
    
    
%     % correlation with gamma
%     px = squeeze(nanmean(nanmean(btf.powspctrm(:,btf.freq>80,:),2)));
%     py = sum(seqevidence{1}(round(linspace(1,size(seqevidence{1},1),length(btf.time))),:),2);
%     pidx = sum(isnan([px py]),2)==0;
%     px = px(pidx);
%     py = py(pidx);
%     [rho,pval] = corr(py,px);
%     
%     % point biserial
%     px = squeeze(nanmean(nanmean(btf.powspctrm(:,btf.freq>80,:),2)));
%     py = otimes';
%     pidx = sum(isnan([px py]),2)==0;
%     px = px(pidx);
%     py = py(pidx);
%     [rho,h,pval] = pointbiserial(py,px);
    
    oidx = find(otimes);
    trialinfo = [];
    for o = 1:length(oidx)
        t0 = findMin(btf.time(oidx(o))-.1,btf.time);
        t1 = findMin(btf.time(oidx(o))+.15,btf.time);
        trialinfo = [trialinfo; t0 t1];
    end
    trialinfo(:,3) = trialinfo(:,2) - trialinfo(:,1);
    trialinfo = trialinfo(trialinfo(:,3)==median(trialinfo(:,3)),:);
    
    if ~isempty(trialinfo)
        replay_tf = [];
        for o = 1:size(trialinfo,1)
            replay_tf(o,:,:) = squeeze(nanmean(tf.powspctrm(:,:,trialinfo(o,1):trialinfo(o,2)))); % average over channels
        end
        if size(trialinfo,1)>1
            replay_tf = replay_tf(~isnan(sum(squeeze(mean(replay_tf,2)),2)),:,:); % remove replay events that go outside the bounds of time-frequency time window
        end
        rTF{trl} = replay_tf;
        
        x = linspace(-.1,.15,size(replay_tf,3));
        replay_btf = replay_tf - squeeze(nanmean(nanmean(replay_tf(:,:,x>=(-.1) & x<=(-.05)),3)));
        rbTF{trl} = replay_btf;
        
        %{
            figure
            imagesc(squeeze(mean(replay_btf)),'XData',linspace(-100,150,size(replay_btf,3)),'YData',tf.freq);
            colormap(colours(100,'inferno'))
            view(180,90);
            set(gca,'xdir','reverse')
            hold on
            plot([0 0],btf.freq([1 end]),'w:');
        %}
    end
end

save(fullfile(directories.dir_save,[subject '_tf.mat']),'TF','bTF','rTF','rbTF','freq','timex','-v7.3');

% Plot
if makePlots
    
    % --- get maximum trial length (in time-frequency samples)
    dur = nan(1,nTrls);
    for trl = 1:nTrls
        dur(trl) = size(TF{trl},2);
    end

    % --- overlap time
    x = tf.freq;
    avtf = nan(nTrls,length(x),max(dur));
    for trl = 1:nTrls
        avtf(trl,:,1:size(TF{trl},2)) = TF{trl};
    end

    figure
    imagesc(squeeze(nanmean(avtf)),'XData',tf.time,'YData',tf.freq);
    colormap(colours(100,'inferno'))
    view(180,90);
    set(gca,'xdir','reverse')
    hold on
    plot([0 0],tf.freq([1 end]),'w:');
    drawnow

    % --- average over time into frequency bands
    freqbands = [-Inf 4
        4 8
        8 12
        12 30
        30 80
        80 Inf];
    fdata = [];
    for c = 1:2
        idx = behav.Choice==c;
        for f = 1:size(freqbands,1)
            tmp = avtf(idx,x>=freqbands(f,1) & x<freqbands(f,2),:);
            fdata(c,f) = nanmean(tmp(:));
        end
    end

    figure
    bar(fdata')
    set(gca,'xticklabels',{'Delta','Theta','Alpha','Beta','Low Gamma','High Gamma'})
    legend({'Approach','Avoid'})
end

thisT = behav(:,ismember(behav.Properties.VariableNames,...
    {'Subject','Practice','Block','Trial','ExpTrial','P','Choice','nV_1','nV_2','RT','Acc','EV'}));
thisT.Theta = squeeze(nanmean(nanmean(avtf(:,x>=4 & x<8,:),3),2)); % average across time
thisT.Alpha = squeeze(nanmean(nanmean(avtf(:,x>=8 & x<12,:),3),2)); % average across time
thisT.Beta = squeeze(nanmean(nanmean(avtf(:,x>=12 & x<30,:),3),2)); % average across time
thisT.LowGamma = squeeze(nanmean(nanmean(avtf(:,x>=30 & x<80,:),3),2)); % average across time
thisT.HighGamma = squeeze(nanmean(nanmean(avtf(:,x>=80,:),3),2)); % average across time

save(fullfile(directories.dir_save,[subject '_tf-table.mat']),'thisT');

if makePlots
    % --- show replay averages
    avreplay = [];
    for trl = 1:nTrls
        if size(rbTF{trl},1)==1
            avreplay(trl,:,:) = squeeze(rbTF{trl});
        else
            avreplay(trl,:,:) = squeeze(nanmean(rbTF{trl})); % average across replay events
        end
    end

    figure
    imagesc(squeeze(nanmean(avreplay)),'XData',linspace(-100,150,size(avreplay,3)),'YData',tf.freq);
    colormap(colours(100,'inferno'))
    view(180,90);
    set(gca,'xdir','reverse')
    hold on
    plot([0 0],tf.freq([1 end]),'w:');
    drawnow
end

end