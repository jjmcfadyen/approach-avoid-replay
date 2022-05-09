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

    [thiscv,thislambda] = get_bestLambdas(subjects{s},trainTimes,thisnull);
    CV(s,:,:) = thiscv;
    lambdas(s,:) = thislambda;
    
end

%% Epoch the replay onsets for each subject (SPM format)

lagrange = 20:10:90; % lags at which to look for onsets (in ms)

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

dir_save = 'D:\2020_RiskyReplay\data\meg\replay\epochs_nullthresh_bc_paths';

epochtype = 'short'; % 'short' = -100 to 150ms, 'long' = -1000 to 1000 ms
switch epochtype
    case 'short'
        twin = [-.3 .3];
    case 'long' % NOTE: THIS PRODUCES A VERY LARGE AMOUNT OF DATA (~300-400 GB)
        twin = [-1 1];
end

for s = 1:N
   
    disp(['Getting replay onsets for ' subjects{s} '...'])
    
    if ~exist(fullfile(dir_save,subjects{s}))
        mkdir(fullfile(dir_save,subjects{s}))
    end
    
    % Get classifier
    thislambda = lambdas(s,trainTimes==optimised_times(s));
    tmp = load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(optimised_times(s)) '_n' num2str(thisnull) '.mat']));
    classifier = tmp.classifier;
    classifier.betas = squeeze(classifier.betas(:,:,thislambda));
    classifier.intercepts = classifier.intercepts(:,thislambda);
    
    % Get task data (100 Hz)
    tmp = load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subjects{s} '_task_100Hz.mat'])); % loads 'merged' variable
    merged = tmp.merged;
    if ~isfield(merged,'fsample')
        merged.fsample = 100;
    end
    
    % Identify replay onsets
    load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
    behav = behav.task;

    params = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    onsets = get_replayOnsets(params,classifier,lagrange/10,'nullperm',behav); % lag range in samples
    
    % Make epochs for the replay onsets per trial
    filelist = params.block;
    fileinclude = ones(length(filelist),1);
    for run = 1:size(params,1)

        D = spm_eeg_load(fullfile(dir_meg,'5_ICA_ds-600Hz',subjects{s},['ICA_ds-600Hz_' subjects{s} '_task_r' num2str(params.block(run)) '.mat']));

        blockonsets = onsets(onsets.Run==params.block(run),:);

        if ~isempty(blockonsets)

            nOnsets = size(blockonsets,1);

            % Epoch replay events
            S = [];
            S.D = D;
            S.bc = 0;
            S.trl = round([(blockonsets.Time - abs(twin(1)))/(1/D.fsample),... % onset, minus 100 ms baseline
                     (blockonsets.Time + abs(twin(2)))/(1/D.fsample),... % offset (onset + 150 ms)
                     repmat(twin(1)*D.fsample,size(blockonsets,1),1)]); % trial shift to accomodate baseline

%             S.conditionlabels = repmat({'anyreplay'},nOnsets,1);
            S.conditionlabels = cell(nOnsets,1);
            S.conditionlabels(blockonsets.Path==1) = {'rewardingreplay'};
            S.conditionlabels(blockonsets.Path==0) = {'aversivereplay'};

            % Epoch
            epoched = spm_eeg_epochs(S);

            % Baseline correct using -100 to -50 ms
            S = [];
            S.D = epoched;
            S.timewin = [-100 -50];
            S.prefix = '';
            epoched = spm_eeg_bc(S);

            % Move file
            epoched.move(fullfile(dir_save,subjects{s},['replay-epochs_r' num2str(params.block(run)) '_' subjects{s} '.mat']));
        else
            fileinclude(run,:) = 0;
        end
    end  
    filelist = filelist(find(fileinclude));
    
    % Merge
    S = [];
    S.D = cell(1,length(filelist));
    for f = 1:length(filelist)
        S.D{f} = spm_eeg_load(fullfile(dir_save,subjects{s},['replay-epochs_r' num2str(filelist(f)) '_' subjects{s} '.mat']));
    end
    if strcmp(subjects{s},'506559')
        S.D(:,1) = []; % fiducials missing for first block
    end
    replayepochs = spm_eeg_merge(S);

    if ~exist(fullfile(dir_save,'merged'))
        mkdir(fullfile(dir_save,'merged'))
    end
    replayepochs.move(fullfile(dir_save,'merged',['replay-epochs_all_' subjects{s} '.mat']));
    
    if any(contains(chantype(replayepochs),'EEG'))
        replayepochs = spm_eeg_load(fullfile(dir_save,'merged',['replay-epochs_all_' subjects{s} '.mat']));
        replayepochs = chantype(replayepochs,find(contains(chantype(replayepochs),'EEG')),'Other');
        replayepochs.save;
    end
end

%% Do time-frequency analysis in SPM

clear matlabbatch
cc = 0;
for s = 1:N

    cc = cc+1;
    matlabbatch{cc}.spm.meeg.tf.tf.D = {['D:\2020_RiskyReplay\data\meg\replay\epochs_nullthresh_bc\merged\replay-epochs_all_' subjects{s} '.mat']};
    matlabbatch{cc}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{cc}.spm.meeg.tf.tf.frequencies = [1:150];
    matlabbatch{cc}.spm.meeg.tf.tf.timewin = [-Inf Inf];
    matlabbatch{cc}.spm.meeg.tf.tf.method.morlet.ncycles = 5;
    matlabbatch{cc}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{cc}.spm.meeg.tf.tf.method.morlet.subsample = 6;
    matlabbatch{cc}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{cc}.spm.meeg.tf.tf.prefix = '';

    cc = cc+1;
    matlabbatch{cc}.spm.meeg.averaging.average.D = {['D:\2020_RiskyReplay\data\meg\replay\epochs_nullthresh_bc\merged\tf_replay-epochs_all_' subjects{s} '.mat']};
    matlabbatch{cc}.spm.meeg.averaging.average.userobust.standard = false;
    matlabbatch{cc}.spm.meeg.averaging.average.plv = false;
    matlabbatch{cc}.spm.meeg.averaging.average.prefix = 'm';
end

spm_jobman('run',matlabbatch);

% Convert to fieldtrip
TF = cell(1,N);
bTF = cell(1,N);
logTF = cell(1,N);
logbTF = cell(1,N);
for s = 1:N

    % load data
    D = spm_eeg_load(['D:\2020_RiskyReplay\data\meg\replay\epochs_nullthresh_bc\merged\mtf_replay-epochs_all_' subjects{s} '.mat']);
    data = ftraw(D,1:length(D.chanlabels),1:length(D.time),1:length(D.condlist));

    % get time frequency estimate across trials in this condition (i)
    TF{s} = data;

    % calculate log power
    logTF{s} = TF{s};
    logTF{s}.powspctrm = log(logTF{s}.powspctrm);

    % baseline correct normal and log-transformed
    cfg = [];
    cfg.baseline = [-.1 -.05];
    bTF{s} = ft_freqbaseline(cfg, TF{s,1});

    cfg = [];
    cfg.baseline = [-.1 -.05];
    logbTF{s} = ft_freqbaseline(cfg, logTF{s});

    %{
    cfg = [];
    cfg.layout = 'CTF275.lay';
    cfg.baseline = [-.1 -.05];
    cfg.xlim = [-.1 .15];
    cfg.ylim = [4 150];
    figure
    ft_multiplotTFR(cfg,TF{s})
    sgtitle(['s=' num2str(s)])
    drawnow
    %}

    clear tmp

end

% plot
grandTF = ft_freqgrandaverage([],logbTF{:,1});

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.baseline = [-.1 -.05];  
cfg.xlim = [-.1 .15];
cfg.ylim = [2 50];
cfg.colormap = colours(256,'viridis');
figure
ft_multiplotTFR(cfg,grandTF)

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.baseline = [-.1 -.05];
cfg.xlim = [-.1 .15];
cfg.ylim = [50 150];
cfg.colormap = colours(256,'viridis');
figure
ft_multiplotTFR(cfg,grandTF)


% Stats
nulldata = cell(1,N);
for s = 1:N
    nulldata{s} = logbTF{s};
    nulldata{s}.powspctrm = zeros(size(nulldata{s}.powspctrm));
end

cfg = [];
cfg.method = 'template';
cfg.layout = 'CTF275.lay';
neighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.channel          = {'MEG'};
cfg.frequency        = [50 150];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 100;

cfg_neighb.method    = 'distance';
cfg.neighbours       = neighbours;

cfg.design           = [ones(1,N) ones(1,N)*2];
cfg.ivar             = 1;

cfg.latency = [0 0.1];
cfg.avgoverfreq = 'yes';

stat = ft_freqstatistics(cfg, logbTF{:}, nulldata{:});

figure
imagesc(squeeze(mean(stat.stat))); %imagesc(1-stat{i}.prob);
colormap('hot')
xlabel('Time')
ylabel('Channels')
set(gca,'ytick',1:length(stat.label));
set(gca,'yticklabels',stat.label);
title(['Lowest p = ' num2str(min(stat.prob(:)))])
set(gcf,'position',[440  -141 484 957])
drawnow;


% Cross-frequency coupling
for s = 1:N

    crossfreq = ft_crossfrequencyanalysis(cfg, freq);

end

%% Do beamforming in OSL (send to cluster)

frequencies = [120 150];

conditionType = 'path'; % 'any' or 'path'

switch conditionType
    case 'any'
        conditions = {'anyreplay'};
        foldername = '';
    case 'path'
        conditions = {'rewardingreplay','aversivereplay'};
        foldername = '_paths';
end

% addpath('D:\Toolboxes\spm12')
% spm('defaults','eeg')

% directories on this work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';
dir_epochs = fullfile(dir_meg,'replay',['epochs_nullthresh_bc' foldername],'merged');

cluster_type = 'holly'; % 'holly' or 'myriad'

% directories on cluster
switch cluster_type
    case 'myriad'
        dir_clustermeg = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/meg/';
    case 'holly'
        dir_clustermeg = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/';
end
dir_continuous = [dir_clustermeg '5_ICA_ds-600Hz/'];
dir_epoch = [dir_clustermeg '8_replayepochs-600Hz_nullthresh_bc' foldername '/merged/'];

for s = 1:N
    
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);

    dir_output = [dir_clustermeg,['beamforming/anyreplay_nullthresh_merged_bc' foldername '_' num2str(frequencies(1)) '-' num2str(frequencies(2)) 'Hz/'],subjects{s}];

    oat = [];
    
    rr = 1;

    epochedfile                     = [dir_epoch,'/replay-epochs_all_' subjects{s} '.mat'];
    oat.source_recon.D_epoched{rr}  = epochedfile; % the above file epoched into replay onsets (-100 to 100 ms), no baseline correction

    oat.source_recon.conditions   = conditions; % as determined by conditionType setting
    oat.source_recon.freq_range   = frequencies;
    oat.source_recon.time_range   = [-.1 .15];

    oat.source_recon.method                         = 'beamform';
    oat.source_recon.normalise_method               = 'mean_eig';
    oat.source_recon.gridstep                       = 5;
    oat.source_recon.forward_meg                    = 'Single Shell';
    oat.source_recon.modalities                     = {'MEG'};
    oat.source_recon.report.do_source_variance_maps = 1;

    oat.source_recon.dirname = dir_output;

    oat.first_level.design_matrix_summary             = {};
    oat.first_level.contrast                          = {};
    oat.first_level.contrast_name                     = {};
    switch conditionType
        case 'any'
            oat.first_level.design_matrix_summary{1}  = [1]; % one condition (replay onset for any path)
            oat.first_level.contrast{1}               = [1];
            oat.first_level.contrast_name{1}          = 'replay';
            
        case 'path'
            oat.first_level.design_matrix_summary{1}  = [1 0]; % rewarding replay
            oat.first_level.design_matrix_summary{2}  = [0 1]; % aversive replay
            
            oat.first_level.contrast{1}               = [1 0];
            oat.first_level.contrast_name{1}          = 'rewarding';
            oat.first_level.contrast{2}               = [0 1];
            oat.first_level.contrast_name{2}          = 'aversive';
            oat.first_level.contrast{3}               = [1 -1];
            oat.first_level.contrast_name{3}          = 'differential';
        case 'choice'
            oat.first_level.design_matrix_summary{1}  = [1 0 0 0]; % avoid:    rewarding replay
            oat.first_level.design_matrix_summary{2}  = [0 1 0 0]; % avoid:    aversive replay
            oat.first_level.design_matrix_summary{3}  = [0 0 1 0]; % approach: rewarding replay
            oat.first_level.design_matrix_summary{4}  = [0 0 0 1]; % approach: aversive replay
            
            oat.first_level.contrast{1}               = [1 -1 0 0];
            oat.first_level.contrast_name{1}          = 'diff_avoid';
            oat.first_level.contrast{2}               = [0 0 1 -1];
            oat.first_level.contrast_name{2}          = 'diff_approach';
            oat.first_level.contrast{3}               = [1 1 -1 -1];
            oat.first_level.contrast_name{3}          = 'diff_choice';
            oat.first_level.contrast{4}               = [1 -1 -1 1];
            oat.first_level.contrast_name{4}          = 'interaction';
    end

    oat.first_level.report.first_level_cons_to_do   = 1;
%     oat.first_level.cope_type                       = 'acope'; % to address sign ambiguity

    oat.first_level.time_range                      = [-.1 .1];
    oat.first_level.baseline_timespan               = [-.1 -.05];
    oat.first_level.name                            = ['replay'];

    oat.subject_level.subjects_to_do = 1;
    oat.subject_level.session_index_list = {1:rr};

    oat.to_do = [1 1 1 0];

    % save
    save(fullfile(dir_batch,['oat_' subjects{s} '.mat']),'oat');
    generate_jobs_beamforming(['oat_' subjects{s} '.mat'],cluster_type);
    
end
disp('done')

%% Do second level in SPM

regressor = ''; % '' for none, 'acc', 'anxiety'
smoothing = 12;

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

replaytype = 'anyreplay_nullthresh_merged_bc_4-8Hz'; 

if contains(replaytype,'paths')
    excludeSubjects = {'263098','680913'};
    C = 3;
else
    excludeSubjects = [];
    C = 1;
end

thesesubjects = setdiff(subjects,excludeSubjects);
thisN = length(thesesubjects);

dir_images = fullfile('D:\2020_RiskyReplay\data\meg\beamforming\',replaytype);
dir_group = fullfile(dir_images,'group');

% --- smooth images
% clear matlabbatch
% cc = 0;
% 
% for c = 1:C
%     cc = cc+1;
%     for s = 1:thisN
%         filename = fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir',['tstat' num2str(c) '_5mm.nii']); 
%         if ~exist(filename) && exist([filename '.gz'])
%             gunzip([filename '.gz']);
%         end
%         matlabbatch{cc}.spm.spatial.smooth.data{s,1} = filename;
%     end
%     matlabbatch{cc}.spm.spatial.smooth.fwhm = repmat(smoothing,1,3);
%     matlabbatch{cc}.spm.spatial.smooth.dtype = 0;
%     matlabbatch{cc}.spm.spatial.smooth.im = 0;
%     matlabbatch{cc}.spm.spatial.smooth.prefix = ['s' num2str(smoothing)];
% end
% 
% spm_jobman('run',matlabbatch);

if ~isempty(regressor)

    regvec = nan(thisN,1);
    switch regressor
        case 'acc'
            for s = 1:thisN
                load(fullfile(dir_behav,thesesubjects{s},[thesesubjects{s} '_parsedBehav.mat']))
                behav = behav.task;
                regvec(s,1) = mean(behav.Acc(behav.Forced==0 & behav.RT<=30));
            end
        case 'anxiety'
            tmp = readtable('D:\2020_RiskyReplay\results\questionnaire_results.csv');
            for s = 1:thisN
                idx = find(tmp.Subject==str2double(thesesubjects{s}));
                regvec(s,1) = mean([tmp.score_IUS(idx) tmp.score_worry(idx)]);
            end
    end
    regvec = regvec - mean(regvec);
end

% GLM on each time sample (for visualisation)
if contains(replaytype,'paths')
    C = 2; % 1 = rewarding vs aversive (paired sample t-test), 2 = differential (one-sample t-test)
else
    excludeSubjects = [];
    C = 1;
end

timepoints = -0.05:0.01:0.15;
for tp = 1:length(timepoints)
    for c = 1:C
        
        if ~(~isempty(regressor) && c==2)
            if c==1
                thisdir = [dir_group '_onesample_' regressor '_c' num2str(c) '_s' num2str(smoothing) '_t' num2str(timepoints(tp)*1000)];
            elseif c==2
                thisdir = [dir_group '_pairedsample_' regressor '_c' num2str(c) '_s' num2str(smoothing) '_t' num2str(timepoints(tp)*1000)];
            end
            if ~exist(thisdir)
                mkdir(thisdir)
            end
    
            clear matlabbatch
            cc = 0;
    
            cc = cc+1;
            matlabbatch{cc}.spm.stats.factorial_design.dir = {thisdir};
            if c==1
                if contains(replaytype,'paths')
                    thiscon = 3; % 'differential' replay contrast
                else
                    thiscon = 1; % 'any' replay contrast
                end
                % (one sample)
                for s = 1:thisN
                    subjectdir = fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir');
                    times = importdata(fullfile(subjectdir,'times'));
                    matlabbatch{cc}.spm.stats.factorial_design.des.t1.scans{s,1} = fullfile(subjectdir,['s' num2str(smoothing) 'tstat' num2str(thiscon) '_5mm.nii,' num2str(findMin(timepoints(tp),times))]);
                end
            elseif c==2
                % (paired samples)
                for s = 1:thisN
                    subjectdir = fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir');
                    times = importdata(fullfile(subjectdir,'times'));
                    matlabbatch{cc}.spm.stats.factorial_design.des.pt.pair(s).scans = {
                        fullfile(subjectdir,['s' num2str(smoothing) 'tstat1_5mm.nii,' num2str(findMin(timepoints(tp),times))])
                        fullfile(subjectdir,['s' num2str(smoothing) 'tstat2_5mm.nii,' num2str(findMin(timepoints(tp),times))])
                        };
                end
            end
            matlabbatch{cc}.spm.stats.factorial_design.des.pt.gmsca = 0;
            matlabbatch{cc}.spm.stats.factorial_design.des.pt.ancova = 0;
            
            if ~isempty(regressor)
                if c==1
                    matlabbatch{cc}.spm.stats.factorial_design.cov.c = regvec;
                elseif c==2
                    matlabbatch{cc}.spm.stats.factorial_design.cov.c = squash([regvec regvec]);
                end
                matlabbatch{cc}.spm.stats.factorial_design.cov.cname = regressor;
                matlabbatch{cc}.spm.stats.factorial_design.cov.iCFI = 1;
                matlabbatch{cc}.spm.stats.factorial_design.cov.iCC = 1;
            end
    
            matlabbatch{cc}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{cc}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{cc}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{cc}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{cc}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{cc}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{cc}.spm.stats.factorial_design.globalm.glonorm = 1;
    
            cc = cc+1;
            matlabbatch{cc}.spm.stats.fmri_est.spmmat = {fullfile(thisdir,'SPM.mat')};
            matlabbatch{cc}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{cc}.spm.stats.fmri_est.method.Classical = 1;
    
            cc = cc+1;
            matlabbatch{cc}.spm.stats.con.spmmat = {fullfile(thisdir,'SPM.mat')};
            if c==1
                if ~isempty(regressor)
                    zerovec = [0];
                else
                    zerovec = [];
                end
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.name = 'pos';
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.weights = [1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.name = 'neg';
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.weights = [-1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.name = 'any';
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.weights = [1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
                if ~isempty(regressor)
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.name = regressor;
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.weights = [0 1];
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                    matlabbatch{cc}.spm.stats.con.consess{5}.fcon.name = regressor;
                    matlabbatch{cc}.spm.stats.con.consess{5}.fcon.weights = [0 1];
                    matlabbatch{cc}.spm.stats.con.consess{5}.fcon.sessrep = 'none';
                end
            elseif c==2
                if ~isempty(regressor)
                    zerovec = zeros(1,thisN+1);
                else
                    zerovec = zeros(1,thisN);
                end
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.name = 'rewarding';
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.weights = [1 -1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.name = 'aversive';
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.weights = [-1 1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.name = 'diff';
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.weights = [1 -1 zerovec];
                matlabbatch{cc}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
                if ~isempty(regressor)
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.name = regressor;
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 zerovec(1:end-1)];
                    matlabbatch{cc}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                    matlabbatch{cc}.spm.stats.con.consess{5}.tcon.name = regressor;
                    matlabbatch{cc}.spm.stats.con.consess{5}.tcon.weights = [0 0 -1 zerovec(1:end-1)];
                    matlabbatch{cc}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
                    matlabbatch{cc}.spm.stats.con.consess{6}.fcon.name = regressor;
                    matlabbatch{cc}.spm.stats.con.consess{6}.fcon.weights = [0 0 1 zerovec(1:end-1)];
                    matlabbatch{cc}.spm.stats.con.consess{6}.fcon.sessrep = 'none';
                end
            end
    
            % Run job(s)
            spm_jobman('run',matlabbatch);
        end
    end
end

% % --- average across 0 to 100 ms window
% if contains(replaytype,'paths')
%     C = 3; % 1 = rewarding vs aversive (paired sample t-test), 2 = differential (one-sample t-test)
% else
%     excludeSubjects = [];
%     C = 1;
% end
% 
avgtime = [0 0.1]; % in seconds, to average across
% clear matlabbatch
% cc = 0;
% for c = 1:C
%     for s = 1:thisN
%         
%         times = importdata(fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir','times'));
%         frames = find(times>=avgtime(1) & times<=avgtime(end));
%         
%         cc = cc+1;
%         expression = '(';
%         for f = 1:length(frames)
%             matlabbatch{cc}.spm.util.imcalc.input{f,1} = fullfile(dir_images,[thesesubjects{s} '.oat'],...
%                 'session1_replay_dir',['s' num2str(smoothing) 'tstat' num2str(c) '_5mm.nii,' num2str(frames(f))]);
%             expression = [expression,'i' num2str(f)];
%             if f<length(frames)
%                 expression = [expression ' + '];
%             else
%                 expression = [expression ')/' num2str(length(frames))];
%             end
%         end
%         matlabbatch{cc}.spm.util.imcalc.expression = expression;
%         matlabbatch{cc}.spm.util.imcalc.output = ['avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms_s' num2str(smoothing) 'tstat' num2str(c) '_5mm'];
%         matlabbatch{cc}.spm.util.imcalc.outdir = {fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir')};
%         matlabbatch{cc}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%         matlabbatch{cc}.spm.util.imcalc.options.dmtx = 0;
%         matlabbatch{cc}.spm.util.imcalc.options.mask = 0;
%         matlabbatch{cc}.spm.util.imcalc.options.interp = 1;
%         matlabbatch{cc}.spm.util.imcalc.options.dtype = 4;
%     end
% end
% spm_jobman('run',matlabbatch);

% GLM on average window   
if contains(replaytype,'paths')
    C = 2; % one sample, two sample
else
    excludeSubjects = [];
    C = 1;
end

for c = 1:C

    if ~(~isempty(regressor) && c==2)
        if c==1
            thisdir = [dir_group '_onesample_c1_' regressor '_s' num2str(smoothing) '_avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms'];
        elseif c==2
            thisdir = [dir_group '_twosample_c2_' regressor '_s' num2str(smoothing) '_avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms'];
        end
        if ~exist(thisdir)
            mkdir(thisdir)
        end
    
        clear matlabbatch
        cc = 0;
        cc = cc+1;
        matlabbatch{cc}.spm.stats.factorial_design.dir = {thisdir};
        if c==1
            if contains(replaytype,'paths')
                thiscon = 3; % 'differential' replay contrast
            else
                thiscon = 1; % 'any' replay contrast
            end
            % (one sample)
            for s = 1:thisN
                subjectdir = fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir');
                matlabbatch{cc}.spm.stats.factorial_design.des.t1.scans{s,1} = fullfile(subjectdir,['avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms_s' num2str(smoothing) 'tstat' num2str(thiscon) '_5mm.nii,1']);
            end
        elseif c==2
            % (paired samples)
            for s = 1:thisN
                subjectdir = fullfile(dir_images,[thesesubjects{s} '.oat'],'session1_replay_dir');
                matlabbatch{cc}.spm.stats.factorial_design.des.pt.pair(s).scans = {
                    fullfile(subjectdir,['avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms_s' num2str(smoothing) 'tstat1_5mm.nii,1'])
                    fullfile(subjectdir,['avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms_s' num2str(smoothing) 'tstat2_5mm.nii,1'])
                    };
            end
        end
        matlabbatch{cc}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{cc}.spm.stats.factorial_design.des.pt.ancova = 0;
        
        if ~isempty(regressor)
            if c==1
                matlabbatch{cc}.spm.stats.factorial_design.cov.c = regvec;
            elseif c==2
                matlabbatch{cc}.spm.stats.factorial_design.cov.c = squash([regvec regvec]);
            end
            matlabbatch{cc}.spm.stats.factorial_design.cov.cname = regressor;
            matlabbatch{cc}.spm.stats.factorial_design.cov.iCFI = 1;
            matlabbatch{cc}.spm.stats.factorial_design.cov.iCC = 1;
        end
        
        matlabbatch{cc}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{cc}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{cc}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{cc}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{cc}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{cc}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{cc}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        cc = cc+1;
        matlabbatch{cc}.spm.stats.fmri_est.spmmat = {fullfile(thisdir,'SPM.mat')};
        matlabbatch{cc}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{cc}.spm.stats.fmri_est.method.Classical = 1;
        
        cc = cc+1;
        matlabbatch{cc}.spm.stats.con.spmmat = {fullfile(thisdir,'SPM.mat')};
        if c==1
            if ~isempty(regressor)
                zerovec = [0];
            else
                zerovec = [];
            end
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.name = 'pos';
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.weights = [1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.name = 'neg';
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.weights = [-1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.name = 'any';
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.weights = [1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
            if ~isempty(regressor)
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.name = regressor;
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.weights = [0 1];
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.name = regressor;
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.weights = [0 1];
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.sessrep = 'none';
            end
        elseif c==2
            if ~isempty(regressor)
                zerovec = zeros(1,thisN+1);
            else
                zerovec = zeros(1,thisN);
            end
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.name = 'rewarding';
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.weights = [1 -1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.name = 'aversive';
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.weights = [-1 1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.name = 'diff';
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.weights = [1 -1 zerovec];
            matlabbatch{cc}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
            if ~isempty(regressor)
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.name = regressor;
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.weights = [0 0 zerovec(1:end-1) 1];
                matlabbatch{cc}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.name = regressor;
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.weights = [0 0 zerovec(1:end-1) 1];
                matlabbatch{cc}.spm.stats.con.consess{5}.fcon.sessrep = 'none';
            end
        end
        
        % Run job(s)
        spm_jobman('run',matlabbatch);
    end
end

%% Plot MNI coordinate over time

% --> Open this the SPM results in the GUI and right clice and select 'Extract table as data structure' - this gives you TabDat
% e.g., D:\2020_RiskyReplay\data\meg\beamforming\anyreplay_4-8Hz\group_onesample\SPM.mat

coord = nan(size(TabDat.dat,1),3);
for i = 1:size(TabDat.dat,1)
    coord(i,:) = TabDat.dat{i,end}';
end

timepoints = -.05:0.01:0.1;
m = nan(1,length(timepoints));
for tp = 1:length(timepoints)
    cd(fullfile('D:\2020_RiskyReplay\data\meg\beamforming\',replaytype,['group_onesample_c1_s' num2str(smoothing) '_t' num2str(timepoints(tp)*1000)]))
    load('SPM.mat');
    
    for c = 1:size(coord,1)
        fname = SPM.xCon(1).Vspm.fname;
        vcoords = mm2vox(coord(c,:),spm_vol(fname));
        m(c,tp) = SPM.xCon(1).Vspm.private.dat(vcoords(1),vcoords(2),vcoords(3));
    end
end

figure
plot(timepoints,m,'linewidth',1.4); hold on
plot(timepoints([1 end]),[0 0],'k:');

% y = [];
% x = [];
% for s = 1:length(SPM.xY.VY)
%     vcoords = mm2vox(coord,spm_vol(SPM.xY.VY(s).fname));
%     V = SPM.xY.VY(s).private;
%     y(s,:) = V.dat(vcoords(1),vcoords(2),vcoords(3),:);
%     x(s,:) = V.timing.toffset:V.timing.tspace:(V.timing.tspace*size(V.dat,4)+V.timing.toffset-V.timing.tspace);
% end
% x = mean(x);
% 
% m = mean(y);
% sem = std(y)/sqrt(size(y,1));
% upper = m+sem;
% lower = m-sem;
% 
% figure
% patch([x fliplr(x)],[upper fliplr(lower)],'k','facealpha',.2,'edgecolor','none'); hold on
% plot(x,m,'k','linewidth',1.4); hold on
% plot(x([1 end]),[0 0],'k:');

%% Beamforming in DAiSS toolbox (on work PC)

dir_save = 'D:\2020_RiskyReplay\data\meg\beamforming';
dir_epochs = 'D:\2020_RiskyReplay\data\meg\replay\epochs_nullthresh\merged';

freqband = [80 256];
conditions = {
            'aversivereplay_avoid'
            'rewardingreplay_avoid'
            'aversivereplay_approach'
            'rewardingreplay_approach'
            }';

for s = 2:N
   
    dir_subj = fullfile(dir_save,subjects{s});
    if ~exist(dir_subj)
        mkdir(dir_subj)
    end
    
    clear matlabbatch
    cc = 0;
    
    % Head model
    cc = cc+1;
    matlabbatch{cc}.spm.meeg.source.headmodel.D = {fullfile(dir_epochs,subjects{s},['replay-epochs_all_' subjects{s} '.mat'])};
    matlabbatch{cc}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{cc}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'FIL_CTF_L';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'FIL_CTF_R';
    matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{cc}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{cc}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    
    % Prepare
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.data.dir = {dir_subj};
    matlabbatch{cc}.spm.tools.beamforming.data.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
    matlabbatch{cc}.spm.tools.beamforming.data.val = 1;
    matlabbatch{cc}.spm.tools.beamforming.data.gradsource = 'inv';
    matlabbatch{cc}.spm.tools.beamforming.data.space = 'MNI-aligned';
    matlabbatch{cc}.spm.tools.beamforming.data.overwrite = 0;
    
    % Define sources
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{cc}.spm.tools.beamforming.sources.reduce_rank = [2 3];
    matlabbatch{cc}.spm.tools.beamforming.sources.keep3d = 1;
    matlabbatch{cc}.spm.tools.beamforming.sources.plugin.grid.resolution = 5;
    matlabbatch{cc}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';
    matlabbatch{cc}.spm.tools.beamforming.sources.plugin.grid.constrain = 'iskull';
    matlabbatch{cc}.spm.tools.beamforming.sources.visualise = 1;
    
    % Covariance features
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.features.BF = {fullfile(dir_subj,'BF.mat')};
    matlabbatch{cc}.spm.tools.beamforming.features.whatconditions.condlabel = conditions;
    matlabbatch{cc}.spm.tools.beamforming.features.woi = [-Inf Inf];
    matlabbatch{cc}.spm.tools.beamforming.features.modality = {'MEG'};
    matlabbatch{cc}.spm.tools.beamforming.features.fuse = 'no';
    matlabbatch{cc}.spm.tools.beamforming.features.plugin.csd.foi = freqband;
    matlabbatch{cc}.spm.tools.beamforming.features.plugin.csd.taper = 'dpss';
    matlabbatch{cc}.spm.tools.beamforming.features.plugin.csd.keepreal = 0;
    matlabbatch{cc}.spm.tools.beamforming.features.plugin.csd.hanning = 1;
    matlabbatch{cc}.spm.tools.beamforming.features.regularisation.minkatrunc.reduce = 1;
    matlabbatch{cc}.spm.tools.beamforming.features.bootstrap = false;
    
    % Inverse solution
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{cc-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{cc}.spm.tools.beamforming.inverse.plugin.dics.fixedori = 'yes';
    
    % Output
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{cc-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.reference.power = 1;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.powmethod = 'lambda1';
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.whatconditions.condlabel = conditions;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.sametrials = false;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.woi = [-Inf Inf];
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.contrast = 1;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.logpower = false;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.foi = freqband;
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.taper = 'dpss';
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.result = 'bycondition'; % 'singleimage' or 'bycondition' or 'bytrial'
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.scale = 'yes';
    matlabbatch{cc}.spm.tools.beamforming.output.plugin.image_dics.modality = 'MEG';
    
    % Write
    cc = cc+1;
    matlabbatch{cc}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{cc-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{cc}.spm.tools.beamforming.write.plugin.nifti.normalise = 'all'; % 'all' to normalise across all images, 'separate' to normalise each image separately
    matlabbatch{cc}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';
    
    % Run
    spm_jobman('run',matlabbatch);
    
end

% Do stats
dir_group = 'D:\2020_RiskyReplay\data\meg\beamforming\group';
if ~exist(dir_group)
    mkdir(dir_group)
end

factors = {'choice','pathtype'};
conditions = { 'aversivereplay_avoid'
    'rewardingreplay_avoid'
    'aversivereplay_approach'
    'rewardingreplay_approach'}; % I think this is the order they appear in but I'm not sure... (this is what it is in D.inv{1}.inverse.trials')
levels = [2 2;
          2 1
          1 2
          1 1];

clear matlabbatch
cc = 0;

cc = cc+1;
matlabbatch{cc}.spm.stats.factorial_design.dir = {dir_group};
for f = 1:length(factors)
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).name = factors{f};
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).levels = 2;
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).dept = 0;
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).variance = 1;
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).gmsca = 0;
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).ancova = 0;
end
for l = 1:size(levels,1)
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.icell(l).levels = levels(l,:);
    for s = 1:N
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.icell(l).scans{s,1} = fullfile(dir_save,subjects{s},...
            ['dics_pow_cond_' conditions{l} '_replay-epochs_all_' subjects{s} '.nii']);
    end
end
matlabbatch{cc}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{cc}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{cc}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{cc}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{cc}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{cc}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{cc}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{cc}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{cc}.spm.stats.factorial_design.globalm.glonorm = 1;

cc = cc+1;
matlabbatch{cc}.spm.stats.fmri_est.spmmat = {fullfile(dir_group,'SPM.mat')};
matlabbatch{cc}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{cc}.spm.stats.fmri_est.method.Classical = 1;

% Run job(s)
spm_jobman('run',matlabbatch);

%% Get average time between onset events

avlag = [];
for s = 1:N
    
    tmp = replay_onsets(contains(replay_onsets.Subject,subjects{s}),:);
    blockonsets = [tmp.Path tmp.Onset_time];
    
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
