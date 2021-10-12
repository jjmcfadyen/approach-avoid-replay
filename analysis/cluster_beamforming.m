function cluster_beamforming(filename)
% Script to preprocess the data
% (needs to be run on UNIX based system using OSL toolbox)
% r = run number (double)

%% Directories

dir_osl = '/data/holly-host/jmcfadyen/osl-core-master-UNIX';
addpath(genpath(dir_osl));
osl_startup;

%% Beamforming (all events)

% Load OAT variable
load(filename);

if length(oat.first_level.contrast) ~= 4
    error('Incorrect OAT file - not the right number of contrasts...')
end

% Coregister first
for r = 1:length(oat.source_recon.D_continuous)
    for e = 1:2
        
        S = [];
        if e==1
            S.D = oat.source_recon.D_continuous{r};
        elseif e==2
            S.D = oat.source_recon.D_epoched{r};
        end
        
        D = spm_eeg_load(S.D);
        if any(contains(D.chantype,'EEG'))
            D = chantype(D,find(contains(D.chantype,'EEG')),'Other');
            D.save;
        end
        
        S.mri = []; % empty means to use template
        S.useheadshape = 0;
        S.use_rhino = 0;
        S.forward_meg = 'Single Shell';
        S.fid.label.nasion = 'nas';
        S.fid.label.lpa = 'lpa';
        S.fid.label.rpa = 'rpa';
        D = osl_headmodel(S);
        D.save;
    end
end

% Run source reconstruction first level
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

oat.to_do = [0 1 0 0];
nRuns = length(oat.source_recon.D_epoched);
includeRun = ones(1,nRuns);
for r = 1:nRuns
    oat.first_level.sessions_to_do = r;
    try
        oat = osl_run_oat(oat);
    catch
        includeRun(r) = 0;
    end
end

% Run subject level
oat.to_do = [0 0 1 0];
oat.first_level.sessions_to_do = find(includeRun);
oat.subject_level.session_index_list = {find(includeRun)};
oat = osl_run_oat(oat);

% Save
save(fullfile(oat.source_recon.dirname,'oat.mat'),'oat');

% Write .nii files for subject level result
S                           = [];
S.oat                       = oat;
S.stats_fname               = oat.subject_level.results_fnames{1}; %[oat.first_level.results_fnames{f} '.mat'];
S.first_level_contrasts     = 1:length(oat.first_level.contrast);
S.resamp_gridstep           = oat.source_recon.gridstep;
[statsdir,times,count]      = oat_save_nii_stats(S);

% % Write SPM.mat files
% for f = 1:length(oat.first_level.results_fnames)
%     S                   = [];
%     S.oat               = oat;
%     S.stats_fname       = oat.first_level.results_fnames{f};
%     [D_tstat, D_cope]   = oat_save_spm_stats(S);
% end

% Delete redundant files
filelist = dir(fullfile(oat.source_recon.dirname,'concat*'));
for f = 1:length(filelist)
    delete(fullfile(filelist(f).folder,filelist(f).name));
end

filelist = dir(fullfile(oat.source_recon.dirname,'efsession*'));
for f = 1:length(filelist)
    delete(fullfile(filelist(f).folder,filelist(f).name));
end

filelist = dir(fullfile(oat.source_recon.dirname,'session*.mat'));
for f = 1:length(filelist)
    delete(fullfile(filelist(f).folder,filelist(f).name));
end

filelist = dir(fullfile(oat.source_recon.dirname,'TF*'));
for f = 1:length(filelist)
    delete(fullfile(filelist(f).folder,filelist(f).name));
end

end
