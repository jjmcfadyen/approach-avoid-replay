function cluster_preprocess(filename)
% Script to preprocess the data 
% (needs to be run on UNIX based system using OSL toolbox)
% r = run number (double)

%% Directories

dir_data = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/meg/';

dir_osl = '~/Scratch/2021_HrvojeProject/toolboxes/osl-core-master-UNIX/';
addpath(genpath(dir_osl));
osl_startup;

%% Set parameters

% Extract info from filename
filesplit = strsplit(filename, '_');
subject = filesplit{2};
task = filesplit{3};
run = strsplit(filesplit{4}(2:end),'.mat');
run = str2double(run{1});

% Decide which frequency (or frequencies) to downsample to
Fs = [100 600]; % downsampling rate (in Hz)

%% Filter
    
% % file setup
% thisdir = fullfile(dir_data,'2_cropped',subject);
% dir_opts = fullfile(thisdir,['opts_' task '_r' num2str(run)]);
% 
% opts                    = [];
% opts.datatype           = 'ctf';
% opts.spm_files          = {fullfile(thisdir,filename)};
% opts.dirname            = dir_opts;
% 
% % filter
% opts.highpass.do        = 1;
% opts.highpass.cutoff    = 0.5;
% opts.mains.do           = 1;
% 
% % switch off other settings
% opts.downsample.do      = 0;
% opts.bad_segments.do    = 0;
% opts.africa.do          = 0;
% opts.outliers.do        = 0;
% opts.coreg.do           = 0;
% opts.epoch.do           = 0;
% 
% % run
% opts = osl_run_opt(opts); % creates new folder, same as dir_output but with '.opt' appended
% 
% disp(['=================================================================='])
% disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files{1} ' ==='])
% disp(['=================================================================='])
% 
% % movefile
% dir_output = fullfile(dir_data,'3_filtered',subject);
% if ~exist(dir_output)
%     mkdir(dir_output)
% end
% 
% D = spm_eeg_load(opts.results.spm_files{1});
% D.move(fullfile(dir_output,['filtered_' subject '_' task '_r' num2str(run)]));
% 
% rmdir(dir_opts,'s')

%% Downsample & detect bad segments

for ds = 1%:length(Fs)
    
%     % file setup
%     thisdir = fullfile(dir_data,'3_filtered',subject);
%     dir_opts = fullfile(thisdir,['opts_' task '_r' num2str(run)]);
% 
%     opts                    = [];
%     opts.datatype           = 'ctf';
%     opts.spm_files          = {fullfile(thisdir,['filtered_' subject '_' task '_r' num2str(run)])};
%     opts.dirname            = dir_opts;
% 
%     % downsampling
%     opts.downsample.do      = 1;
%     opts.downsample.freq    = Fs(ds);
% 
%     % bad segments
%     opts.bad_segments.do    = 1;
% 
%     % switch off
%     opts.africa.do          = 0;
%     opts.epoch.do           = 0;
%     opts.outliers.do        = 0;
%     opts.coreg.do           = 0;
%     opts.highpass.do        = 0;
%     opts.mains.do           = 0;
% 
%     % run
%     opts = osl_run_opt(opts);
% 
%     disp(['=================================================================='])
%     disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files{1} ' ==='])
%     disp(['=================================================================='])
% 
%     % movefile
%     dir_output = fullfile(dir_data,['4_downsampled_' num2str(Fs(ds)) 'Hz'],subject);
%     if ~exist(dir_output)
%         mkdir(dir_output)
%     end
% 
%     D = spm_eeg_load(opts.results.spm_files{1});
%     D.move(fullfile(dir_output,['ds-' num2str(Fs(ds)) 'Hz_' subject '_' task '_r' num2str(run)]));
% 
%     rmdir(dir_opts,'s')

    %% AFRICA

    % reset random number generator seed for ICA
    rng(str2double(subject));

    % file setup
    thisdir = fullfile(dir_data,['4_downsampled_' num2str(Fs(ds)) 'Hz'],subject);
    dir_opts = fullfile(thisdir,['opts_' task '_r' num2str(run)]);

    opts                    = [];
    opts.datatype           = 'ctf';
    opts.spm_files          = {fullfile(thisdir,['ds-' num2str(Fs(ds)) 'Hz_' subject '_' task '_r' num2str(run)])};
    opts.dirname            = dir_opts;

    % ICA
    opts.africa.do                              = 1;
    opts.africa.todo.ica                        = 1;
    opts.africa.todo.ident                      = 1;
    opts.africa.todo.remove                     = 1;
    opts.africa.ident.max_num_artefact_comps    = 20;
    opts.africa.ident.artefact_chans            = {'EOG'};
    opts.africa.precompute_topos                = false;
    opts.africa.ident.do_kurt                   = true;
    opts.africa.ident.mains_kurt_thresh         = 0.5;
    opts.africa.ident.do_cardiac                = true;
    opts.africa.ident.do_plots                  = false;
    opts.africa.ident.do_mains                  = false;

    % switch off
    opts.downsample.do                          = 0;
    opts.highpass.do                            = 0;
    opts.bad_segments.do                        = 0;
    opts.epoch.do                               = 0;
    opts.outliers.do                            = 0;
    opts.coreg.do                               = 0;

    % run
    opts = osl_run_opt(opts);

    disp(['=================================================================='])
    disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files{1} ' ==='])
    disp(['=================================================================='])

    % movefile
    dir_output = fullfile(dir_data,['5_ICA_ds-' num2str(Fs(ds)) 'Hz'],subject);
    if ~exist(dir_output)
        mkdir(dir_output)
    end

    D = spm_eeg_load(opts.results.spm_files{1});
    D.move(fullfile(dir_output,['ICA_ds-' num2str(Fs(ds)) 'Hz_' subject '_' task '_r' num2str(run)]));

    rmdir(dir_opts,'s')

    %% Epoch
% 
%     if ~grep(scans.task(r),'rest')
% 
%         % file setup
%         opts.spm_files = opts.results.spm_files;
%         opts.datatype  = 'ctf';
%         opts.dirname   = dir_output;
% 
%         % epoch
%         opts.epoch.do    = 1;
%         opts.outliers.do = 1;
% 
%         % switch off
%         opts.coreg.do           = 0;
%         opts.bad_segments.do    = 0;
%         opts.downsample.do      = 0;
%         opts.africa.do          = 0;
%         opts.africa.todo.ica    = 0;
%         opts.africa.todo.ident  = 0;
%         opts.africa.todo.remove = 0;
% 
%         % loop for each event type in this run
%         opts.epoch.trialdef = struct();
%         if grep(scans.task(r),'size')
% 
%             % Size task:
%             % --- 1. Picture 1 (either the feature on its own, or
%             %                   the combined feature stimulus)
%             % --- 2. Picture 2 (either the feature on its own, or
%             %                   the combined feature stimulus)
%             % --- 3. Size estimate
%             % --- 4. Feedback
%             %
%             % The 'phase' column in T indicates whether the feature
%             % or stimulus was presented first.
%             %
%             % Value of trials is the IMAGE ID (e.g., "Darius3.png")
%             % that participants were estimating the size of (for both
%             % the probe and the stimulus)
% 
%             thisorder = unique(rT.phase);
%             if strcmp(thisorder{1},'probe-stimulus')
%                 epoch_labels = {'probe-first','stimulus-second'};
%             elseif strcmp(thisorder{1},'stimulus-probe')
%                 epoch_labels = {'stimulus-first','probe-second'};
%             end
% 
%             epoch_values = unique(rT.probe); % probes (1 to 4)
%             stimtype = unique(rT.stimImage); % name of stimulus (e.g. "Darius")
% 
%             idx = find(~contains(extractfield(D.events,'type'),'artefact_OSL'));
%             thesetrials = extractfield(D.events,'value');
%             thesetrials = thesetrials(idx);
%             if iscell(thesetrials)
%                 thesetrials = cell2mat(thesetrials);
%             end
%             thesetrials = unique(thesetrials);
% 
%             theseprobes = T.probe(thesetrials);
% 
%             cc = 0;
%             for i = 1:length(epoch_labels) % for stimulus and probe...
% 
%                 eventtype = strsplit(epoch_labels{i},'-');
%                 eventtype = eventtype{1};
% 
%                 for j = 1:length(epoch_values) % for probes 1 to 4...
%                     cc = cc + 1;
%                     opts.epoch.trialdef(cc).conditionlabel  = [epoch_labels{i} '-' stimtype{1} '-' epoch_values{j}]; % e.g. 'probe-first-3'
%                     opts.epoch.trialdef(cc).eventtype        = eventtype; % e.g. 'probe' (matching ev.type)
%                     opts.epoch.trialdef(cc).eventvalue       = thesetrials(strfindcell(theseprobes,epoch_values{j})); % matching trials (matching ev.value)
%                     opts.epoch.time_range                   = [-0.1 0.4]; 
%                 end
%             end
% 
%             opts = osl_run_opt(opts);
% 
%             disp(['=================================================================='])
%             disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files_epoched{1} ' ==='])
%             disp(['=================================================================='])
% 
%             % copy file to results folder
%             D = spm_eeg_load(opts.results.spm_files_epoched{1});
%             D = D.copy(fullfile(dir_results,'3_epochs',[filename '_epoched.mat']));
% 
%             disp(['=================================================================='])
%             disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'3_epochs',[filename '_epoched.mat']) ' ==='])
%             disp(['=================================================================='])
% 
%         elseif grep(scans.task(r),'fnc')
% 
%             % Function learning task:
%             % --- 1. Stimulus (combined feature stimulus - 0.4 seconds)
%             % --- 2. Deliberation (2 seconds)
%             % --- 3. Feedback (LEARNING ONLY - absent in TEST)
%             %
%             % We have to do the deliberation & stimuli separately as
%             % they have different time ranges
% 
%             ev_type = extractfield(D.events,'type');
%             ev_value = extractfield(D.events,'value'); % trial number
%             thesetrials = ev_value(contains(ev_type,'stimulus') & ~contains(ev_type,'artefact_OSL'));
%             if iscell(thesetrials)
%                 thesetrials = cell2mat(thesetrials);
%             end
% 
%             for c = 1:2
% 
%                 opts.epoch.trialdef = struct(); % reset
% 
%                 if c == 1
%                     opts.epoch.trialdef(1).conditionlabel   = 'stimulus'; 
%                     opts.epoch.trialdef(1).eventtype        = 'stimulus'; 
%                     opts.epoch.trialdef(1).eventvalue       = thesetrials; % matching trials (matching ev.value)
%                     opts.epoch.time_range                   = [-0.1 0.4]; 
%                 elseif c == 2
%                     opts.epoch.trialdef(1).conditionlabel   = 'deliberation'; 
%                     opts.epoch.trialdef(1).eventtype        = 'deliberation'; 
%                     opts.epoch.trialdef(1).eventvalue       = thesetrials; % matching trials (matching ev.value)
%                     opts.epoch.time_range                   = [-0.1 2]; 
%                 end
% 
%                 opts = osl_run_opt(opts);
% 
%                 disp(['=================================================================='])
%                 disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files_epoched{1} ' ==='])
%                 disp(['=================================================================='])
% 
%                 % copy file to results folder
%                 D = spm_eeg_load(opts.results.spm_files_epoched{1});
% 
%                 disp(['=================================================================='])
%                 if c == 1
%                     D = D.copy(fullfile(dir_results,'3_epochs',[filename '_stimuli_epoched.mat']));
%                     disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'3_epochs',[filename '_stimuli_epoched.mat']) ' ==='])
%                 elseif c == 2
%                     D = D.copy(fullfile(dir_results,'3_epochs',[filename '_deliberation_epoched.mat']));
%                     disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'3_epochs',[filename '_deliberation_epoched.mat']) ' ==='])
%                 end
%                 disp(['=================================================================='])
% 
%             end
%         end
%     end
% 
%     %% Organise files
% 
%     % Delete output folder and its contents
%     rmdir(dir_output,'s');
% end
% 
% %% Merge runs per task type
% 
% if r == nRuns+1
%     
%     INFO = cell(1,2); % one cell per task (size, function)
%     for t = 1:length(tasktypes)
% 
%         INFO{t} = struct;
%         S = struct;
% 
%         if t == 1
%             E = 1; % only one epoch length for stimulus & probe
%         elseif t == 2
%             E = 2; % stimulus epochs are 500ms, deliberation epochs are 2000ms
%         end
% 
%         for e = 1:E
% 
%             if t==1
%                 filenames = dir(fullfile(dir_results,'3_epochs',['*' tasktypes{t} '*_epoched.mat']));
%             elseif t==2 && e==1
%                 filenames = dir(fullfile(dir_results,'3_epochs',['*' tasktypes{t} '*stimuli_epoched.mat']));
%             elseif t==2 && e==2
%                 filenames = dir(fullfile(dir_results,'3_epochs',['*' tasktypes{t} '*deliberation_epoched.mat']));
%             end
% 
%             S.D = cell(length(filenames),1);
%             for f = 1:length(filenames)
% 
%                 % load data
%                 S.D{f,1} = spm_eeg_load(fullfile(filenames(f).folder,filenames(f).name));
% 
%                 % save information
%                 thistask = strsplit(filenames(f).name,'task-');
%                 thistask = strsplit(thistask{2}, '_run');
%                 thistask = thistask{1};
% 
%                 if e==1 && f==1
%                     cc = 1;
%                 else
%                     cc = cc+1;
%                 end
% 
%                 INFO{t}(cc).task = thistask;
%                 INFO{t}(cc).fname = S.D{f,1}.fname;
%                 INFO{t}(cc).badchannels = S.D{f,1}.badchannels;
%                 INFO{t}(cc).badtrials = S.D{f,1}.badtrials;
%                 INFO{t}(cc).conditions = unique(S.D{f,1}.conditions);
% 
%                 Tidx = nan(1,length(S.D{f,1}.events)); % get row number in T that each epoch corresponds to
%                 for i = 1:length(S.D{f,1}.events)
%                     idx = find(~contains(extractfield(S.D{f,1}.events{i},'type'),'artefact_OSL'),1,'first');
%                     Tidx(i) = S.D{f,1}.events{i}(idx).value;
%                 end
% 
%                 INFO{t}(cc).Tidx = Tidx;
% 
%             end
% 
%             D = spm_eeg_merge(S);
%             oldfile = D.fullfile;
% 
%             % Copy & delete
%             if t==1
%                 filename = [subject '_size.mat'];
%             elseif t==2
%                 if e==1
%                     filename = [subject '_fnc_stimuli.mat'];
%                 elseif e==2
%                     filename = [subject '_fnc_deliberation.mat'];
%                 end
%             end
%             D = D.copy(fullfile(dir_results,'4_merged',filename));
% 
%             D = spm_eeg_load(oldfile);
%             D.delete;
%             clear D;
% 
%         end
%     end
% 
%     save(fullfile(dir_results,'4_merged',[subject '_INFO.mat']),'INFO');
% 
%     %% Visualise conditions
% 
%     for t = 1:length(tasktypes)
% 
%         if t==1
%             E=1;
%         elseif t==2
%             E=2;
%         end
% 
%         for e = 1:E
% 
%             if t==1
%                 filename = [subject '_size'];
%             elseif t==2
%                 if e==1
%                     filename = [subject '_fnc_stimuli'];
%                 elseif e==2
%                     filename = [subject '_fnc_deliberation'];
%                 end
%             end
% 
%             S = struct;
%             S.D = spm_eeg_load(fullfile(dir_results,'4_merged',[filename '.mat']));
%             S.prefix = 'avg_';
%             D = spm_eeg_average(S);
%             
%             conditionlabels = D.conditions';
%             data = ftraw(D);
%             
%             save(fullfile(dir_results,'4_merged',[filename '_ft.mat']),'data');
% 
% %             restoredefaultpath
% %             addpath(dir_ft);
% %             ft_defaults;
% % 
% %             % Look at each condition
% %             timelock = cell(1,length(conditionlabels));
% %             for c = 1:length(conditionlabels)
% %                 cfg = [];
% %                 cfg.trials = find(contains(conditionlabels,conditionlabels{c}));
% %                 timelock{c} = ft_timelockanalysis(cfg,data);
% %             end
% %             
% %             if t==1
% %                 J = 2;
% %             elseif t==2
% %                 J = 1;
% %             end
% %             
% %             for j = 1:J
% %             
% %                 if t==1
% %                     if j==1
% %                         idx = contains(conditionlabels,'probe');
% %                     elseif j==2
% %                         idx = contains(conditionlabels,'stimulus');
% %                     end
% %                     cmap = colours(sum(idx),'viridis');
% %                     thisdata = timelock(idx);
% %                 elseif t==2
% %                     if j==1
% %                         cmap = [255 165 51]/255; % orange (stimuli)
% %                     elseif j==2
% %                         cmap = [137 111 253]/255; % purple (deliberation)
% %                     end
% %                     thisdata = timelock;
% %                 end
% %             
% %                 figure
% %                 cfg = [];
% %                 cfg.layout = 'CTF275.lay';
% %                 cfg.linecolor = cmap;
% %                 cfg.linewidth = 1.6;
% %                 
% %                 ft_multiplotER(cfg,thisdata{:});
% %                 
% %                 if t==1
% %                     grand = ft_timelockgrandaverage([],timelock{:});
% %                 elseif t==2
% %                     grand = timelock{1};
% %                 end
% % 
% %                 figure
% %                 cfg = [];
% %                 cfg.layout = 'CTF275.lay';
% %                 cfg.xlim = [.1 .1];
% %                 ft_topoplotER(cfg,grand);
% %             end
%         end
%     end
% 
% end
%     
end