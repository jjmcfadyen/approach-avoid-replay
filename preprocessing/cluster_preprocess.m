function cluster_preprocess(filename)
% Script to preprocess the data 
% (needs to be run on UNIX based system using OSL toolbox)
% r = run number (double)

%% Directories

dir_data = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/meg/';
dir_behav = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/data/behav/';

dir_osl = '~/Scratch/2020_RiskyReplay/toolboxes/osl-core-master-UNIX/';
addpath(genpath(dir_osl));
osl_startup;

%% Set parameters

% Extract info from filename
[subject, task, run] = split_filename(filename);

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

for ds = 2%1:length(Fs)
    
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
    if ~exist(dir_opts)
        mkdir(dir_opts);
    end

    opts                    = [];
    opts.datatype           = 'ctf';
    opts.spm_files          = {fullfile(thisdir,['ds-' num2str(Fs(ds)) 'Hz_' subject '_' task '_r' num2str(run)])};
    opts.dirname            = dir_opts;
    
    % (check that EOG channels have been defined)
    D = spm_eeg_load(opts.spm_files{1});
    if sum(contains(D.chantype,'EOG')) ~= 3
        error([filename ' does not have 3 EOG channels'])
    end

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
    
end
end
