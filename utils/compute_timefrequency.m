function thisT = compute_timefrequency(subject,optimised_time,directories,waveletwidth)
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

% Get maximum trial length
trial_res = 0.01; % time step size
trialdur = nan(nTrls,1);
for trl = 1:nTrls
    trialdur(trl) = length(data.time{trl}(1):trial_res:data.time{trl}(end));
end

% Do time-frequency for each trial
TF = [];
replay_onsets = [];
for trl = 1:nTrls % DO NOT DO THIS IN PARALLEL - the time-frequency computation uses LOADS of CPU (50% usage)
    
    % Get time-frequency across entire trial (600 Hz data)
    disp(['||| TRIAL ' num2str(trl) ' of ' num2str(nTrls) ' |||'])
    
    cfg = [];
    cfg.trials = trl;
    thistrial = ft_selectdata(cfg,data);
    thistrial.trialinfo = data.trialinfo(trl,:);
    
    cfg = [];
    cfg.channel    = 'MEG';
    cfg.method     = 'wavelet';
    cfg.width      = waveletwidth;
    cfg.output     = 'pow';
    cfg.foi        = 1:1:150;
    cfg.toi        = thistrial.time{1}(1):trial_res:thistrial.time{1}(end);
    cfg.keeptrials = 'no';
    cfg.pad        = 'nextpow2';
    
    tf = ft_freqanalysis(cfg, thistrial);
    if isempty(TF)
        TF = tf;
        if trialdur(trl) ~= max(trialdur)
            TF.time = data.time{find(trialdur==max(trialdur),1,'first')}(1):trial_res:data.time{find(trialdur==max(trialdur),1,'first')}(end);
            TF.powspctrm = nan(nTrls,size(TF.powspctrm,2),max(trialdur));
        end
        TF = rmfield(TF,'cfg'); % to save on file size
    end
    TF.powspctrm(trl,:,1:length(tf.time)) = squeeze(mean(tf.powspctrm)); % average across channels (otherwise file sizes are MASSIVE)
    
    % Get replay onsets from the 100 Hz file
    cfg = [];
    cfg.trials = trl;
    thislow = ft_selectdata(cfg,lowdata); % 100 hz data for replay
    thislow.trialinfo = thislow.trialinfo(1,:);
    
    [onsets, ~] = get_replayOnsets(thislow,classifier,lagrange);

    % Save
    replay_onsets = [replay_onsets; onsets];
    
end

save(fullfile(directories.dir_save,[subject '_w' num2str(waveletwidth) '_tf.mat']),'TF','replay_onsets','behav','-v7.3');

end