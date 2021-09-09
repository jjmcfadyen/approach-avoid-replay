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

dir_save = 'D:\2020_RiskyReplay\data\meg\replay\epochs';

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
        
        blockonsets = blockonsets(blockonsets.include==1,:); % remove forced-choice & catch trials
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
        
    end  
    
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
    
    replayepochs.move(fullfile(dir_save,subjects{s},['replay-epochs_all_' subjects{s} '.mat']));

    if strcmp(subjects{s},'707132')
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

% % Back up original templates
% mkdir(fullfile(dir_canonical,'orig'));
% filelist = dir(fullfile(dir_canonical));
% for f = 3:length(filelist)
%     copyfile(fullfile(filelist(f).folder,filelist(f).name),...
%     fullfile(dir_canonical,'orig',filelist(f).name));
% end

for s = 1:N
    
    fname = ['replay-epochs_all_' subjects{s} '.mat'];
    
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

        % Generate priors
        if i==2
            mesh = gifti(fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
            meshcoords = spm_eeg_inv_transform_points(mesh.mat,mesh.vertices);
            % match the vertices with cortical vs hipp vs amg
        end
        
        dir_save = fullfile(dir_epochs,subjects{s},['inv' num2str(i)]);
        if ~exist(dir_save)
            mkdir(dir_save)
        end
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
    
%     % Plot results (compare inversions)
%     figure
%     for i = 1:2
% 
%         % average across trial types
%         dir_save = fullfile(dir_epochs,subjects{s},['inv' num2str(i)]);
%         filelist = dir(fullfile(dir_save,['replay-epochs_all_' subjects{s} '_' num2str(i) '_*.gii']));
%         for f = 1:length(filelist)
%             
%             % extract data from stupid SPM file (from spm_mesh_render.m)
%             M = gifti(fullfile(filelist(f).folder,filelist(f).name)); 
%             MM = M;
%             M  = gifti(fullfile(dir_save,MM.private.metadata(1).value));
%             M.cdata = full(MM.cdata);
%             
%             if f==1
%                 mesh = M;
%             else
%                 mesh.faces(:,:,f) = M.faces; 
%                 mesh.vertices(:,:,f) = M.vertices;
%                 mesh.cdata(:,2) = M.cdata;
%             end
%             
%         end
%         mesh.faces = squeeze(mean(mesh.faces,3));
%         mesh.vertices = squeeze(mean(mesh.vertices,3));
%         mesh.cdata = squeeze(mean(mesh.cdata,2));
%         mesh = export(mesh,'patch');
%         
%         % get ready to plot
%         P = struct('vertices',mesh.vertices, 'faces',mesh.faces);
%         curv = spm_mesh_curvature(P) > 0;
%         curv = 0.5 * repmat(curv,1,3) + 0.3 * repmat(~curv,1,3);
%         T = mesh.facevertexcdata';
%         
% %         col = colours(256,'viridis');
%         cm = [255 255 255;
%             154, 33, 255]/255;
%         col = colours(256,cm);
% 
%         cmap = zeros(size(T,2),3);
%         clim = [false NaN NaN];
%         mi = clim(2);
%         ma = clim(3);
%         if size(col,1)>3 && size(col,1) ~= size(T,1)
%             if size(T,1) == 1
%                 if ~clim(1)
%                     mi = min(T(:));
%                     ma = max(T(:));
%                 end
%                 cmap = squeeze(ind2rgb(floor(((T(:)-mi)/(ma-mi))*size(col,1)),col));
%             elseif isequal(size(T),[size(curv,1) 3])
%                 cmap = T; T = T';
%             else
%                 if ~clim(1), mi = min(T(:)); ma = max(T(:)); end
%                 for i=1:size(T,1)
%                     cmap = cmap + squeeze(ind2rgb(floor(((T(i,:)-mi)/(ma-mi))*size(col,1)),col));
%                 end
%             end
%         else
%             if ~clim(1), ma = max(T(:)); end
%             for i=1:size(T,1)
%                 cmap = cmap + T(i,:)'/ma * col(i,:);
%             end
%         end
%         cmap = repmat(~any(T,1),3,1)' .* curv + repmat(any(T,1),3,1)' .* cmap;
%        
%         % make plot
% %         cmap(:,4) = ones(size(cmap,1),1);
% %         patch('vertices',P.vertices,'faces',P.faces,'edgecolor','none',...
% %             'facecolor',[.9 .9 .9],'facealpha',.1,...
% %             'SpecularStrength',.25,'AmbientStrength',0.3,'DiffuseStrength',0.7,'SpecularExponent',10,...
% %             'Clipping','off'); hold on
% %         
% %         cmap(:,4) = normalise(T,.9999,1);
% %         patch('vertices',P.vertices,'faces',P.faces,'edgecolor','none',...
% %                 'facecolor','interp','FaceVertexCData',cmap(:,1:3),...
% %                 'facealpha','interp','FaceVertexAlphaData',cmap(:,4),...
% %                 'SpecularStrength',1,'AmbientStrength',0.3,'DiffuseStrength',0.7,'SpecularExponent',10,...
% %                 'Clipping','off'); hold on
% 
%         if i==1
%             cc = 0;
%         elseif i==2
%             cc = 3;
%         end
%         for tv = 1:3
%             
%             subplot(2,3,tv+cc)
%             
%             patch('vertices',P.vertices,'faces',P.faces,'edgecolor','none',...
%                     'facecolor','interp','FaceVertexCData',cmap(:,1:3),...
%                     'SpecularStrength',1,'AmbientStrength',0.3,'DiffuseStrength',0.7,'SpecularExponent',10,...
%                     'Clipping','off'); hold on
%             if i==2
%                
%                hipp_mesh = gifti('D:\2020_RiskyReplay\data\mri\bothHipp.trans.gii');
%                patch('vertices',hipp_mesh.vertices,'faces',hipp_mesh.faces,'facecolor','none','edgecolor','k','edgealpha',.15); hold on
%                
%                amg_mesh = gifti('D:\2020_RiskyReplay\data\mri\bothAmg.trans.gii');
%                patch('vertices',amg_mesh.vertices,'faces',amg_mesh.faces,'facecolor','none','edgecolor','k','edgealpha',.15); hold on
%                
%                hipp_idx = ismember(mesh.vertices(:,1),hipp_mesh.vertices(:,1)) & ...
%                    ismember(mesh.vertices(:,2),hipp_mesh.vertices(:,2)) & ...
%                    ismember(mesh.vertices(:,3),hipp_mesh.vertices(:,3));
%                amg_idx = ismember(mesh.vertices(:,1),amg_mesh.vertices(:,1)) & ...
%                    ismember(mesh.vertices(:,2),amg_mesh.vertices(:,2)) & ...
%                    ismember(mesh.vertices(:,3),amg_mesh.vertices(:,3));
%                title(['Hipp = ' num2str(round(max(T(hipp_idx)) / max(T) * 100,2)) '% max, ',...
%                    ' Amg = ' num2str(round(max(T(amg_idx)) / max(T) * 100,2)) '% max'])
%                
%             end
%                 
%             axis equal
%             if tv==1
%                 view(180,-90)
%             elseif tv==2
%                 view(90,-15)
%             elseif tv==3
%                 view(-90,0)
%             end
%             camlight('headlight')
%             lighting gouraud
%         end
%     end
end

% Group Inversion
for i = 1%:2

    S = [];
    for s = 1:N
        S = [S; fullfile(dir_epochs,subjects{s},['inv' num2str(i)],['replay-epochs_all_' subjects{s} '.mat'])];
    end

    spm_eeg_inv_group(S);

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

for s = 1:N
    
    subjectparams = parameters(contains(parameters.schar,subjects{s}) & contains(parameters.task,'task'),:);
    nRuns = size(subjectparams,1);

    dir_output = [dir_clustermeg,'beamforming/',subjects{s}];

    oat = [];

    for r = 1:nRuns

        thisrun = subjectparams.block(r);
        continuousfile = [dir_continuous,subjects{s},...
            '/ICA_ds-600Hz_' subjects{s} '_task_r' num2str(subjectparams.block(r)) '.mat'];
        epochedfile = [dir_epoch,subjects{s},...
            '/replay-epochs_r' num2str(subjectparams.block(r)) '_' subjects{s} '.mat'];
        
        D = spm_eeg_load(fullfile(dir_meg,'replay','epochs',subjects{s},...
            ['replay-epochs_r' num2str(subjectparams.block(r)) '_' subjects{s} '.mat']));
        
        oat.source_recon.D_continuous{r} = continuousfile; % continuous data for the block - preprocessed (low-pass filter, downsampled to 100Hz, ICA artefact removal)
        oat.source_recon.D_epoched{r}    = epochedfile; % the above file epoched into replay onsets (-100 to 100 ms), no baseline correction
        
        epochinfo = D.history;
        epochinfo = epochinfo(contains(extractfield(epochinfo,'fun'),'spm_eeg_epochs')).args;
        
        oat.source_recon.epochinfo{r}                  = [];
        oat.source_recon.epochinfo{r}.trl              = epochinfo.trl;
        oat.source_recon.epochinfo{r}.conditionlabels  = epochinfo.conditionlabels;
        oat.source_recon.epochinfo{r}.padding          = epochinfo.eventpadding;
    end

    oat.source_recon.conditions   = {'anyreplay_avoid','anyreplay_approach'}; % do all replay types to start with
    oat.source_recon.freq_range   = [1 150];
    oat.source_recon.time_range   = [-.1 .15];

    oat.source_recon.method                         = 'beamform';
    oat.source_recon.normalise_method               = 'mean_eig';
    oat.source_recon.gridstep                       = 5;
    oat.source_recon.forward_meg                    = 'Single Shell';
    oat.source_recon.modalities                     = {'MEG'};
    oat.source_recon.report.do_source_variance_maps = 1;

    oat.source_recon.dirname = dir_output;

    design_matrix_summary                   = {};
    design_matrix_summary{1}                = [1 1]; % one condition (replay onset for any path)
    oat.first_level.design_matrix_summary   = design_matrix_summary;

    oat.first_level.contrast                        = {};
    oat.first_level.contrast{1}                     = [1];
    oat.first_level.contrast_name                   = {};
    oat.first_level.contrast_name{1}                = 'replay';
    oat.first_level.report.first_level_cons_to_do   = 1;
    oat.first_level.cope_type                       = 'acope'; % to address sign ambiguity

    oat.first_level.time_range                      = [-.1 .1];
    oat.first_level.baseline_timespan               = [-.1 -.05];
    oat.first_level.name                            = ['wholebrain_first_level'];

    oat.subject_level.subjects_to_do = 1;
    oat.subject_level.session_index_list = {1:nRuns};

    oat.to_do = [1 1 1 0];

    % save
    save(fullfile(dir_batch,['oat_' subjects{s} '.mat']),'oat');
    generate_jobs_beamforming(['oat_' subjects{s} '.mat'],cluster_type);
    
end

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
