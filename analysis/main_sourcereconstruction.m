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

dir_save = 'D:\2020_RiskyReplay\data\meg\replay\epochs_unnormalised';

epochtype = 'short'; % 'short' = -100 to 150ms, 'long' = -1000 to 1000 ms
switch epochtype
    case 'short'
        twin = [-.1 .15];
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
    load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(optimised_times(s)) '_n' num2str(thisnull) '.mat']));
    classifier.betas = squeeze(classifier.betas(:,:,thislambda));
    classifier.intercepts = classifier.intercepts(:,thislambda);
    
    % Get task data (100 Hz)
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subjects{s} '_task_100Hz.mat'])); % loads 'merged' variable
    if ~isfield(merged,'fsample')
        merged.fsample = 100;
    end
    
    % Identify replay onsets
    onsets = get_replayOnsets(merged,classifier,lagrange);
    
    % Remove onsets that occur during baseline
    onsets = onsets(onsets.Onset_time>0,:);
    
    % Match with behavioural data
    load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
    behav = behav.task;
    
    % Make epochs for the replay onsets per trial
    filelist = parameters.block(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'));
    fileinclude = ones(length(filelist),1);
    for f = 1:length(filelist)
    
        % Load 600 Hz continuous data for each block
        D = spm_eeg_load(fullfile(dir_meg,'5_ICA_ds-600Hz',subjects{s},['ICA_ds-600Hz_' subjects{s} '_task_r' num2str(filelist(f)) '.mat']));

        % Get events from this continuous run
        events = D.events;
        
        etypes = extractfield(events,'type');
        etime = extractfield(events,'time');
        
        evalues = nan(length(events),1);
        for i = 1:length(evalues)
            if ~ischar(events(i).value)
                evalues(i) = events(i).value;
            end
        end
        etypes = etypes(~isnan(evalues));
        etime = etime(~isnan(evalues));
        evalues = evalues(~isnan(evalues));
        
        nEvents = length(evalues);

        % Align each trial event with replay onsets
        trials = unique(evalues);
        blockonsets = [];
        for trl = 1:length(trials)
            
            % block number * 100 + trial number (e.g., block 2, trial 6 = 206)
            thisval = sprintf('%04d',trials(trl));
            thisblock = str2double(thisval(1:2));
            thistrial = str2double(thisval(3:end));
            
            thispractice=0;
            if thisblock==0
                thispractice=1;
                thisblock=1;
            end
            
            idx = behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial;
            
            decisiononset = etime(evalues==trials(trl) & contains(etypes,'decision')');
            if isempty(decisiononset)
                error('Trial onset not found in continuous data')
            end
            
            trialonsets = onsets(onsets.Practice==thispractice & onsets.Block==thisblock & onsets.Trial==thistrial,:);
            n = size(trialonsets,1);
            
            if n>0
                trialonsets.ctime = trialonsets.Onset_time + decisiononset;
                trialonsets.csample = nan(n,1);
                for i = 1:n
                    trialonsets.csample(i) = D.indsample(trialonsets.ctime(i)); 
                end
                if any(trialonsets.csample==0)
                    error('Sample out of range')
                end

                trialonsets.choice = repmat(behav.Choice(idx),n,1);
                trialonsets.include = repmat(behav.Forced(idx)==0 & ((behav.nV_1(idx)>0 & behav.nV_2(idx)<0) | (behav.nV_1(idx)<0 & behav.nV_2(idx)>0)),n,1);

                if behav.nV_1(idx) > behav.nV_2(idx)
                    rewpath = 1;
                    avpath = 2;
                else
                    rewpath = 2;
                    avpath = 1;
                end
                trialonsets.type(trialonsets.Path==0) = {'any'};
                trialonsets.type(trialonsets.Path==1 & rewpath==1) = {'rewarding'};
                trialonsets.type(trialonsets.Path==2 & rewpath==2) = {'rewarding'};
                trialonsets.type(trialonsets.Path==1 & avpath==1) = {'aversive'};
                trialonsets.type(trialonsets.Path==2 & avpath==2) = {'aversive'};

                blockonsets = [blockonsets; trialonsets];
            end
        end
        
        if ~isempty(blockonsets)
            blockonsets = blockonsets(blockonsets.include==1,:); % remove forced-choice & catch trials
        end
        
        if ~isempty(blockonsets)
            nOnsets = size(blockonsets,1);

            % Epoch replay events
            S = [];
            S.D = D;
            S.bc = 0;
            S.trl = [blockonsets.csample - abs(twin(1))*D.fsample,... % onset, minus 100 ms baseline
                     blockonsets.csample + abs(twin(2))*D.fsample,... % offset (onset + 150 ms)
                     repmat(twin(1)*D.fsample,size(blockonsets,1),1)]; % trial shift to accomodate baseline

            S.conditionlabels = cell(nOnsets,1);
            for i = 1:nOnsets
                if blockonsets.choice(i)==1
                    choicelabel = 'approach';
                else
                    choicelabel = 'avoid';
                end
                pathlabel = [blockonsets.type{i} 'replay_'];
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

            % Move file
            epoched.move(fullfile(dir_save,subjects{s},['replay-epochs_r' num2str(filelist(f)) '_' subjects{s} '.mat']));
        else
            fileinclude(f,:) = 0;
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
        replayepochs = chantype(replayepochs,find(contains(chantype(replayepochs),'EEG')),'Other');
        replayepochs.save;
    end
    
end

%% MSP

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

dir_epochs = 'D:\2020_RiskyReplay\data\meg\replay\epochs\merged';
dir_canonical = 'D:\Toolboxes\spm12\canonical';
dir_newmesh = 'D:\2020_RiskyReplay\data\mri';

% Parameters
meshres = 2; % 1 = coarse, 2 = normal, 3 = fine

if meshres==1
    meshunits = 5124;
elseif meshres==2
    meshunits = 8196;
elseif meshres==3
    meshunits = 20484;
end

for i = 1%:2 % 1 = template cortical mesh, 2 = freesurfer colin cortical mesh
    
    % Overwrite template with either the template or a new mesh
    if i==1
        copyfile(fullfile(dir_canonical,'orig',['cortex_' num2str(meshunits) '.surf.gii']),...
            fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
        copyfile(fullfile(dir_canonical,'orig','single_subj_T1.nii'),...
            fullfile(dir_canonical,'single_subj_T1.nii'));
    elseif i==2
        copyfile(fullfile(dir_canonical,'new','canonicalHippAmg.gii'),...
            fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
        copyfile(fullfile(dir_canonical,'new','orig.nii'),...
            fullfile(dir_canonical,'single_subj_T1.nii'));
    end
    
%     % Generate priors
%     if i==2
%         mesh = gifti(fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
%         meshcoords = spm_eeg_inv_transform_points(mesh.mat,mesh.vertices);
%         % match the vertices with cortical vs hipp vs amg
%     end
    
    dir_save = fullfile(dir_epochs,['inv' num2str(i)]);
    if ~exist(dir_save)
        mkdir(dir_save)
    end
    
    for s = 1:N
        
        fname = ['replay-epochs_all_' subjects{s} '.mat'];
        D = spm_eeg_load(fullfile(dir_epochs,subjects{s},fname));
        filename = fullfile(dir_save,fname);
        if ~exist(filename)
            copy(D,filename);
        end
        
        matlabbatch = {};
        
        % Head model specification
        cc = 1;
        matlabbatch{cc}.spm.meeg.source.headmodel.D = {filename};
        matlabbatch{cc}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{cc}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
        matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshres = 2;
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
        matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
        matlabbatch{cc}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
        matlabbatch{cc}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
        
        % Source inversion
        cc = cc+1;
        matlabbatch{cc}.spm.meeg.source.invert.D = {filename};
        matlabbatch{cc}.spm.meeg.source.invert.val = 1;
        matlabbatch{cc}.spm.meeg.source.invert.whatconditions.condlabel = {
            'aversivereplay_avoid'
            'rewardingreplay_avoid'
            'aversivereplay_approach'
            'rewardingreplay_approach'
            }';
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.woi = [-100 150];
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
        matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
        matlabbatch{cc}.spm.meeg.source.invert.modality = {'MEG'};
        
        % Get results
        cc = cc+1;
        matlabbatch{cc}.spm.meeg.source.results.D = {filename};
        matlabbatch{cc}.spm.meeg.source.results.val = 1;
        matlabbatch{cc}.spm.meeg.source.results.woi = [-100 150];
        matlabbatch{cc}.spm.meeg.source.results.foi = [0 0];
        matlabbatch{cc}.spm.meeg.source.results.ctype = 'evoked';
        matlabbatch{cc}.spm.meeg.source.results.space = 1;
        matlabbatch{cc}.spm.meeg.source.results.format = 'mesh';
        matlabbatch{cc}.spm.meeg.source.results.smoothing = 8;
        
        % Run job(s)
        spm_jobman('run',matlabbatch);
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

%% Do beamforming in OSL (send to cluster)

% addpath('D:\Toolboxes\spm12')
% spm('defaults','eeg')

% directories on this work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

cluster_type = 'holly'; % 'holly' or 'myriad'

% directories on cluster
switch cluster_type
    case 'myriad'
        dir_clustermeg = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/meg/';
    case 'holly'
        dir_clustermeg = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/';
end
dir_continuous = [dir_clustermeg '5_ICA_ds-600Hz/'];
dir_epoch = [dir_clustermeg '8_replayepochs-600Hz/'];

conditionType = 'choice'; % 'any' or 'path'
switch conditionType
    case 'any'
        conditions = {'anyreplay_approach','anyreplay_avoid'};
    case {'path','choice'}
        conditions = {'rewardingreplay_avoid','rewardingreplay_approach','aversivereplay_avoid','aversivereplay_approach'};
end

for s = 1:N
    
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);

    dir_output = [dir_clustermeg,'beamforming/',subjects{s}];

    oat = [];

    rr = 0;
    for r = 1:nRuns

        includeRun = true;
        
        thisrun = subjectparams.block(r);
        continuousfile = [dir_continuous,subjects{s},...
            '/ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat'];
        epochedfile = [dir_epoch,subjects{s},...
            '/replay-epochs_r' num2str(subjectparams.block(r)) '_' subjects{s} '.mat'];
        
        Dcont = spm_eeg_load(fullfile('D:\2020_RiskyReplay\data\meg\5_ICA_ds-600Hz',subjects{s},...
            ['ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat']));
        
        try
            
            D = spm_eeg_load(fullfile(dir_meg,'replay','epochs',subjects{s},...
                ['replay-epochs_r' num2str(subjectparams.block(r)) '_' subjects{s} '.mat']));
        
            epochinfo = D.history;
            epochinfo = epochinfo(contains(extractfield(epochinfo,'fun'),'spm_eeg_epochs')).args;

            % check fiducials are present in both files, and all conditions are available
            if isempty(D.fiducials) || isempty(Dcont.fiducials) || ~all(ismember(conditions,D.condlist))
                includeRun = false;
            end
        catch
            includeRun = false;
        end
        
        if includeRun
            
            rr = rr+1;
            
            oat.source_recon.D_continuous{rr} = continuousfile; % continuous data for the block - preprocessed (low-pass filter, downsampled to 100Hz, ICA artefact removal)
            oat.source_recon.D_epoched{rr}    = epochedfile; % the above file epoched into replay onsets (-100 to 100 ms), no baseline correction

            oat.source_recon.epochinfo{rr}                  = [];
            oat.source_recon.epochinfo{rr}.trl              = epochinfo.trl;
            oat.source_recon.epochinfo{rr}.conditionlabels  = epochinfo.conditionlabels;
            oat.source_recon.epochinfo{rr}.padding          = epochinfo.eventpadding;
        end
    end

    oat.source_recon.conditions   = conditions; % as determined by conditionType setting
    oat.source_recon.freq_range   = [1 150];
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
            oat.first_level.design_matrix_summary{1}  = [1 1]; % one condition (replay onset for any path)
            oat.first_level.contrast{1}               = [1];
            oat.first_level.contrast_name{1}          = 'replay';
            
        case 'path'
            oat.first_level.design_matrix_summary{1}  = [1 1 0 0]; % rewarding replay
            oat.first_level.design_matrix_summary{2}  = [0 0 1 1]; % aversive replay
            
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
    oat.first_level.name                            = ['differential'];

    oat.subject_level.subjects_to_do = 1;
    oat.subject_level.session_index_list = {1:rr};

    oat.to_do = [1 1 1 0];

    % save
    save(fullfile(dir_batch,['oat_' subjects{s} '.mat']),'oat');
    generate_jobs_beamforming(['oat_' subjects{s} '.mat'],cluster_type);
    
end


%% Do second level in SPM

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

replaytype = 'choice_1-150Hz'; % 'anyreplay_4-8Hz'    'differential_1-150Hz'

if contains(replaytype,'differential') || contains(replaytype,'choice')
    firstlevelname = 'differential';
    excludeSubjects = {'263098','680913'};
else
    firstlevelname = 'wholebrain_first_level';
    excludeSubjects = {'945092'}; 
end

thesesubjects = setdiff(subjects,excludeSubjects);
thisN = length(thesesubjects);

dir_images = fullfile('D:\2020_RiskyReplay\data\meg\beamforming\',replaytype);
dir_group = fullfile(dir_images,'group');

% --- smooth images
clear matlabbatch
cc = 0;

for c = 1:4
    cc = cc+1;
    for s = 1:thisN
        filename = fullfile(dir_images,[thesesubjects{s} '.oat'],['subject1_' firstlevelname '_sub_level_dir'],['tstat' num2str(c) '_5mm.nii']); 
        if ~exist(filename) && exist([filename '.gz'])
            gunzip([filename '.gz']);
        end
        matlabbatch{cc}.spm.spatial.smooth.data{s,1} = filename;
    end
    matlabbatch{cc}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{cc}.spm.spatial.smooth.dtype = 0;
    matlabbatch{cc}.spm.spatial.smooth.im = 0;
    matlabbatch{cc}.spm.spatial.smooth.prefix = 's';
end

spm_jobman('run',matlabbatch);

% % --- average across 0 to 100 ms window
% avgtime = [0 0.1]; % in seconds, to average across
% for s = 1:thisN
%     
%     times = importdata(fullfile(dir_images,[thesesubjects{s} '.oat'],'subject1_wholebrain_first_level_sub_level_dir','times'));
%     frames = find(times>=avgtime(1) & times<=avgtime(end));
%     
%     cc = cc+1;
%     expression = '(';
%     for f = 1:length(frames)
%         matlabbatch{cc}.spm.util.imcalc.input{f,1} = fullfile(dir_images,[thesesubjects{s} '.oat'],...
%             'subject1_wholebrain_first_level_sub_level_dir',['ststat1_5mm.nii,' num2str(frames(f))]);
%         expression = [expression,'i' num2str(f)];
%         if f<length(frames)
%             expression = [expression ' + '];
%         else
%             expression = [expression ')/' num2str(length(frames))];
%         end
%     end
%     matlabbatch{cc}.spm.util.imcalc.expression = expression;
%     matlabbatch{cc}.spm.util.imcalc.output = ['avg' num2str(avgtime(1)*1000) '-' num2str(avgtime(end)*1000) 'ms_ststat1_5mm'];
%     matlabbatch{cc}.spm.util.imcalc.outdir = {fullfile(dir_images,[thesesubjects{s} '.oat'],'subject1_wholebrain_first_level_sub_level_dir')};
%     matlabbatch{cc}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{cc}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{cc}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{cc}.spm.util.imcalc.options.interp = 1;
%     matlabbatch{cc}.spm.util.imcalc.options.dtype = 4;
% end

timepoints = -0.1:0.01:0.1;
for tp = 1:length(timepoints)
    for c = 1:4
        
        thisdir = [dir_group '_onesample_c' num2str(c) '_t' num2str(timepoints(tp)*1000)];
        if ~exist(thisdir)
            mkdir(thisdir)
        end

        clear matlabbatch
        cc = 0;

        cc = cc+1;
        matlabbatch{cc}.spm.stats.factorial_design.dir = {thisdir};
        % (one sample)
        for s = 1:thisN
            subjectdir = fullfile(dir_images,[thesesubjects{s} '.oat'],['subject1_' firstlevelname '_sub_level_dir']);
            times = importdata(fullfile(subjectdir,'times'));
            matlabbatch{cc}.spm.stats.factorial_design.des.t1.scans{s,1} = fullfile(subjectdir,['ststat' num2str(c) '_5mm.nii,' num2str(findMin(timepoints(tp),times))]);
        end
        matlabbatch{cc}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{cc}.spm.stats.factorial_design.des.pt.ancova = 0;
        matlabbatch{cc}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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
        matlabbatch{cc}.spm.stats.con.consess{1}.tcon.name = 'pos';
        matlabbatch{cc}.spm.stats.con.consess{1}.tcon.weights = [1];
        matlabbatch{cc}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{cc}.spm.stats.con.consess{2}.tcon.name = 'neg';
        matlabbatch{cc}.spm.stats.con.consess{2}.tcon.weights = [-1];
        matlabbatch{cc}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{cc}.spm.stats.con.consess{3}.fcon.name = 'any';
        matlabbatch{cc}.spm.stats.con.consess{3}.fcon.weights = [1];
        matlabbatch{cc}.spm.stats.con.consess{3}.fcon.sessrep = 'none';

        % Run job(s)
        spm_jobman('run',matlabbatch);
    end
end


% Plot MNI coordinate over time
timepoints = -0.05:0.01:0.1;
% --> Open this the SPM results in the GUI and right clice and select 'Extract table as data structure' - this gives you TabDat
coord = nan(size(TabDat.dat,1),3);
for i = 1:size(TabDat.dat,1)
    coord(i,:) = TabDat.dat{i,end}';
end

m = nan(1,length(timepoints));
for tp = 1:length(timepoints)
    cd(fullfile('D:\2020_RiskyReplay\data\meg\beamforming\',replaytype,['group_onesample_t' num2str(timepoints(tp)*1000)]))
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
dir_epochs = 'D:\2020_RiskyReplay\data\meg\replay\epochs\merged';

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
