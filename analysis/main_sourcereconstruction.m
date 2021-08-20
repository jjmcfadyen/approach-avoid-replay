%% Source reconstruction

clear all
clc

%% Parameters & directories

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

% Classification parameters
trainTimes = 0:10:300; % in ms
nT = length(trainTimes);
thisnull = 1; % proportion of data to replicate as null (zeros)
nStates = 6;

% Optimisation parameters
load(fullfile(dir_replay,'optimised_times.mat'));

%% Get lambda for each classifier

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
        CV(s,t,1) = min(cv(:,best_lambda)); 
        CV(s,t,2) = mean(cv(:,best_lambda)); 
        CV(s,t,3) = max(cv(:,best_lambda)); 
        
    end
end

%% Get replay onsets for each subject

lagrange = 20:10:90; % lags at which to look for onsets (in ms)

replay_onsets = [];
for s = 1:N
   
    disp(['Getting replay onsets for ' subjects{s} '...'])
    
    % Get classifier
    thislambda = lambdas(s,trainTimes==optimised_times(s));
    load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(optimised_times(s)) '_n' num2str(thisnull) '.mat']));
    classifier.betas = squeeze(classifier.betas(:,:,thislambda));
    classifier.intercepts = classifier.intercepts(:,thislambda);
    
    % Get task data
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subjects{s} '_task_100Hz.mat'])); % loads 'merged' variable
    if ~isfield(merged,'fsample')
        merged.fsample = 100;
    end
    
    % Identify replay onsets
    onsets = get_replayOnsets(merged,classifier,lagrange);
    
    % Remove onsets from baseline
    onsets = onsets(onsets.Onset_time>0,:);
    
    % Match with behavioural data
    load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
    behav = behav.task;
    
    % Create table
    onsets.Choice = nan(size(onsets,1),1);
    onsets.Forced = nan(size(onsets,1),1);
    onsets.P = nan(size(onsets,1),1);
    onsets.nV_1 = nan(size(onsets,1),1);
    onsets.nV_2 = nan(size(onsets,1),1);
    onsets.ExpTrial = nan(size(onsets,1),1);
    onsets.EV = nan(size(onsets,1),1);
    onsets.RT = nan(size(onsets,1),1);
    for trl = 1:size(behav,1)
        idx = find(behav.Practice(trl)==onsets.Practice & ...
                   behav.Block(trl)==onsets.Block & ...
                   behav.Trial(trl)==onsets.Trial);
        onsets.Choice(idx)  = behav.Choice(trl);
        onsets.Forced(idx)          = behav.Forced(trl);
        onsets.P(idx)           = behav.P(trl);
        onsets.nV_1(idx)        = behav.nV_1(trl);
        onsets.nV_2(idx)        = behav.nV_2(trl);
        onsets.ExpTrial(idx)    = behav.ExpTrial(trl);
        onsets.EV(idx)          = behav.EV(trl);
        onsets.RT(idx)          = behav.RT(trl);
    end
    
    % Remove onests in post-decision 500 ms window
    onsets(onsets.Onset_time > onsets.RT,:) = [];
    
    % Add to group table
    onsets.Subject = repmat({subjects{s}},size(onsets,1),1);
    replay_onsets = [replay_onsets; onsets];
    
end

%% Epoch the data to the replay onsets

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

epochtype = 'long'; % 'short' = -100 to 150ms, 'long' = -1000 to 1000 ms
switch epochtype
    case 'short'
        twin = [-.1 .15];
    case 'long'
        twin = [-1 1];
end

for s = 1:N
   
    dir_output = fullfile(dir_meg,['8_replayepochs-600Hz'],subjects{s});
    if ~exist(dir_output)
        mkdir(dir_output)
    end
   
    % Load data used to determine replay onsets
    load(fullfile(dir_meg,'7_merged_ds-100Hz',[subjects{s} '_task_100Hz.mat'])); % loads 'merged' variable
    
    % Get behavioural logs for main task
    behav = parse_behav(subjects{s},dir_behav);
    
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);
    
    for r = 1:nRuns
    
        % Get continuous data 
        D = spm_eeg_load(fullfile(dir_meg,'5_ICA_ds-600Hz',subjects{s},...
            ['ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat'])); % loads 'merged' variable
        
        % Get replay onsets for this block
        thisrun = subjectparams.block(r);
        if thisrun==0
            thisblock = 1;
            thispractice = 1;
        else
            thisblock = thisrun;
            thispractice = 0;
        end
        thesereplayonsets = replay_onsets(contains(replay_onsets.Subject,subjects{s}) & replay_onsets.Practice==thispractice & replay_onsets.Block==thisblock,:);
        
        % Convert the replay onsets to continuous time (rather than the time since the start of each trial)
        trls = unique(thesereplayonsets.Trial); % trials in this block
        for i = 1:length(trls)
            idx = thesereplayonsets.Trial==trls(i); % which replay onsets belong to this trial
            sampleonsets = merged.sampleinfo(merged.trialinfo(:,1)==thisrun & merged.trialinfo(:,2)==trls(i),:);
            thesereplayonsets.Onset_sample(idx) = thesereplayonsets.Onset_sample(idx) + sampleonsets(1);
        end
        thesereplayonsets.Onset_time = thesereplayonsets.Onset_sample/100;
        
        if thesereplayonsets.Onset_time(end) > D.time(end)
            error('Replay onsets end AFTER the continuous data!')
        end
        
        % Epoch replay events
        S = [];
        S.D = fullfile(dir_meg,'5_ICA_ds-600Hz',subjects{s},...
            ['ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat']);
        S.bc = 0;
        S.trl = [thesereplayonsets.Onset_sample - abs(twin(1))*D.fsample,... % onset, minus 100 ms baseline
                 thesereplayonsets.Onset_sample + abs(twin(2))*D.fsample,... % offset (onset + 150 ms)
                 repmat(twin(1)*D.fsample,size(thesereplayonsets,1),1)]; % trial shift to accomodate baseline
             
        S.conditionlabels = cell(size(thesereplayonsets,1),1);
        for i = 1:length(S.conditionlabels)
            if thesereplayonsets.Choice(i)==1
                choicelabel = 'approach';
            else
                choicelabel = 'avoid';
            end
            if thesereplayonsets.Path(i)==1 && thesereplayonsets.nV_1(i)>thesereplayonsets.nV_2(i)
                pathlabel = 'rewardingreplay_';
            elseif thesereplayonsets.Path(i)==2 && thesereplayonsets.nV_1(i)<thesereplayonsets.nV_2(i)
                pathlabel = 'rewardingreplay_';
            elseif thesereplayonsets.Path(i)==1 && thesereplayonsets.nV_1(i)<thesereplayonsets.nV_2(i)
                pathlabel = 'aversivereplay_';
            elseif thesereplayonsets.Path(i)==2 && thesereplayonsets.nV_1(i)>thesereplayonsets.nV_2(i)
                pathlabel = 'aversivereplay_';
            end
            S.conditionlabels{i} = [pathlabel choicelabel];
        end
        
        % Epoch
        epoched = spm_eeg_epochs(S);
        
        % Baseline correct using -100 to -50 ms
        S = [];
        S.D = epoched;
        S.timewin = [-100 -50];
        S.prefix = '';
        epoched = spm_eeg_bc(S);
        
        % Move
        epoched.move(fullfile(dir_output,['replay-epochs-' epochtype '_ds-600Hz_' subjects{s} '_task_r' num2str(thisrun) '.mat']));
    
    end
end

%% Do time-frequency analysis in Fieldtrip

epochtype = 'long';

for s = 1:N
   
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);
    
    for r = 1:nRuns
    
        % Get continuous data 
        D = spm_eeg_load(fullfile(dir_meg,'8_replayepochs-600Hz',subjects{s},...
            ['replay-epochs-' epochtype '_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat']));
        
        % Convert to fieldtrip format
        data = ftraw(D);
        
        cfg = [];
        cfg.channel = 'MEG';
        data = ft_selectdata(cfg,data);
        
        % Plot epochs
        cfg = [];
        cfg.layout = 'CTF275.lay';
        cfg.xlim = [-.1 .15];
        figure
        ft_multiplotER(cfg,data);
        
        % Power spectra
        cfg = [];
        cfg.toilim = [-.1 .15];
        shortdata = ft_redefinetrial(cfg,data);
        
        cfg = [];
        cfg.output = 'pow';
        cfg.channel = 'MEG';
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.foi = 20:180;
        PS = ft_freqanalysis(cfg,shortdata);
        
        figure;
        plot(PS.freq,PS.powspctrm); hold on
        plot(PS.freq,mean(PS.powspctrm),'k','linewidth',2);
        
        % Hanning window
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'MEG';
%         cfg.method       = 'mtmconvol';
%         cfg.taper        = 'hanning';
%         cfg.foi          = 2:150;                           % frequencies (in Hz)
%         cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window = 0.1 sec
%         cfg.toi          = -0.1:0.01:0.15;                % time window "slides" from -0.01 to 0.15 sec in steps of 0.01 sec (10 ms)
        cfg.method       = 'wavelet';
        cfg.width        = 7;
        cfg.foi          = 2:150;
        cfg.toi          = -0.1:0.01:0.15;
        TF = ft_freqanalysis(cfg, data);
        
        cfg = [];
        cfg.baseline = [-.1 -.05];
        bTF = ft_freqbaseline(cfg, TF);
        
        cfg = [];
        cfg.layout = 'CTF275.lay';
        cfg.baseline = [-.1 -.05];
        cfg.xlim = [-.1 .15];
        figure
        ft_multiplotTFR(cfg,TF)
        
        figure
        subplot(1,2,1)
        imagesc(flipud(squeeze(mean(TF.powspctrm)))); % average over channels
        subplot(1,2,2)
        imagesc(flipud(squeeze(mean(bTF.powspctrm)))); % average over channels
        
    end
end

%% Do source reconstruction using MSP in SPM

for s = 1:N

    
end

%% Do beamforming in OSL (send to cluster)

% directories on this work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

% directories on cluster
dir_clustermeg = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/meg/';
dir_continuous = [dir_clustermeg '5_ICA_ds-600Hz/'];
dir_epoch = [dir_clustermeg '8_replayepochs-600Hz/'];

for s = 1:N
    
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);
    
    dir_output = [dir_clustermeg,'beamforming/',subjects{s}];
    
    for r = 1:nRuns
        
        thisrun = subjectparams.block(r);
        
        continuousfile = [dir_continuous,subjects{s},...
            '/ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat'];
        epochedfile = [dir_epoch,subjects{s},...
            '/replay-epochs_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat'];
        
        D = spm_eeg_load(fullfile(dir_meg,'8_replayepochs-600Hz',subjects{s},...
            ['replay-epochs_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat']));
        
        oat = [];

        oat.source_recon.D_continuous = {continuousfile}; % continuous data for the block - preprocessed (low-pass filter, downsampled to 100Hz, ICA artefact removal)
        oat.source_recon.D_epoched    = {epochedfile}; % the above file epoched into replay onsets (-100 to 100 ms), no baseline correction
        
        oat.source_recon.conditions   = D.condlist; % do all replay types to start with
        oat.source_recon.freq_range   = [1 150];
        oat.source_recon.time_range   = [-.1 .15];

        oat.source_recon.method                         = 'beamform';
        oat.source_recon.normalise_method               = 'mean_eig';
        oat.source_recon.gridstep                       = 5;
        oat.source_recon.forward_meg                    = 'Single Shell';
        oat.source_recon.modalities                     = {'MEGGRAD'};
        oat.source_recon.report.do_source_variance_maps = 1;

        oat.source_recon.dirname = dir_output;

        design_matrix_summary                   = {};
        design_matrix_summary{1}                = [1]; % one condition (replay onset for any path)
        oat.first_level.design_matrix_summary   = design_matrix_summary;

        oat.first_level.contrast                        = {};
        oat.first_level.contrast{1}                     = [1];
        oat.first_level.contrast_name                   = {};
        oat.first_level.contrast_name{1}                = 'replay';
        oat.first_level.report.first_level_cons_to_do   = 1;

        oat.first_level.time_range                      = [-.1 .1];
        oat.first_level.name                            = ['wholebrain_first_level'];

        oat.to_do = [1 1 0 0];

        % save
        save(fullfile(dir_batch,['oat_' subjects{s} '_r' num2str(thisrun) '.mat']),'oat');
        generate_jobs_beamforming(['oat_' subjects{s} '_r' num2str(thisrun) '.mat']);
        
    end
end

%% Get average time between onset events

avlag = [];
for s = 1:N
    
    tmp = replay_onsets(contains(replay_onsets.Subject,subjects{s}),:);
    theseonsets = [tmp.Path tmp.Onset_time];
    
    % chunk into path replay sections
    thislag = [];
    i = 0;
    while true
       
        i = i+1;

        if i > size(tmp,1)
            break
        end
        
        thispath = tmp.Path(i);
        
        offset = find(tmp.Path((i+1):end) ~= thispath,1,'first')+i;
        if ~isempty(offset) 
            if length(unique(tmp.Trial(i:offset)))>1% don't span trials
                i = offset;
            elseif tmp.Onset_time(offset)-tmp.Onset_time(i) < 0
                break
            else
                thislag = [thislag; tmp.Onset_time(offset)-tmp.Onset_time(i) tmp.Choice(i)];
                i = offset;
            end
        else
            break
        end
    end
    avlag(s,:) = [mean(thislag(thislag(:,2)==1)) mean(thislag(thislag(:,2)==2))];
end

%% Get sequenceness of sequenceness

maxLag = 100; % in samples (100 Hz)
bins = maxLag;
nIterations = 10;

seqseq = cell(1,N);
for s = 1:N
    
    disp(['Getting sequenceness-of-sequenceness for ' subjects{s} '...'])

    % Get classifier
    thislambda = lambdas(s,trainTimes==optimised_times(s));
    load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(optimised_times(s)) '_n' num2str(thisnull) '.mat']));
    classifier.betas = squeeze(classifier.betas(:,:,thislambda));
    classifier.intercepts = classifier.intercepts(:,thislambda);
    
    % Get task data
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subjects{s} '_task_100Hz.mat'])); % loads 'merged' variable
    if ~isfield(merged,'fsample')
        merged.fsample = 100;
    end
    
    % Get evidence for each replay onset, per path
    [~,TM] = get_replayOnsets(merged,classifier,lagrange);
    
    % Do GLM per trial
    thisseq = nan(length(merged.trial),nIterations+1,maxLag);
    parfor trl = 1:length(merged.trial)
       
        xidx = merged.time{trl} >= 0;
        Y = TM{trl}(xidx,:);
        
        nanidx = sum(isnan(Y),2)>0;
        Y = Y(~nanidx,:);
        
        nSamples = size(Y,1);
        nCol = size(Y,2);
        
        % Design matrix - transitions in this permutation
        dM = [0 0 1 1
              0 0 1 1
              1 1 0 0
              1 1 0 0];
        dM = dM(:);

        % add autocorrelation and constant
        dM = [dM squash(eye(nCol)) squash(ones(nCol))];
        
        tmpseq = nan(1,maxLag);
        for it = 1:nIterations+1
            
            disp(['--- getting sequenceness-of-sequenceness for trial ' num2str(trl) ', iteration ' num2str(it) '...'])
            
            if it==1
                thisY = Y;
            else
                thisY = Y(randperm(nSamples),:);
            end
            
            % Create toeplitz matrix
            TP = nan(nSamples,nCol*maxLag);
            cc = linspace(0,size(TP,2),nCol+1);
            warning off
            for st = 1:nCol
                tmp = toeplitz(thisY(:,st),zeros(maxLag+1,1));
                TP(:,cc(st)+1:cc(st+1)) = tmp(:,2:end);
            end
            warning on

            % First-level GLM
            betas = nan(nCol*maxLag, nCol);
            for ilag = 1:bins
                idx = (1:bins:nCol*maxLag) + ilag - 1;
                tmp = pinv([TP(:,idx) ones(size(TP,1),1)])*thisY;
                betas(idx,:) = tmp(1:end-1,:);
            end

            betas = reshape(betas,[maxLag nCol^2]);

            % Second-level GLM
            SEQ = pinv(dM) * betas';

            % Save to variable
            tmpseq(it,:) = SEQ(1,:);
        end
        thisseq(trl,:,:) = tmpseq;
    end
    seqseq{s} = thisseq;
    
    % Plot sequenceness-of-sequenceness
    %{
    figure
    opts = [];
    opts.g = 1;
    ss_plot(seqseq{s},opts);
    %}

end

% Plot all
pd = [];
for s = 1:N
    pd(s,:,:) = squeeze(mean(seqseq{s}));
end

figure
opts = [];
opts.g = 1;
ss_plot(pd,opts);
