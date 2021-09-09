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

% Coregister first
for r = 1:length(oat.source_recon.D_continuous)
    for e = 1:2
        S = [];
        if e==1
            S.D = oat.source_recon.D_continuous{r};
        elseif e==2
            S.D = oat.source_recon.D_epoched{r};
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

% Run source reconstruction
oat = osl_check_oat(oat);
oat = osl_run_oat(oat);

% Save
save(fullfile(oat.source_recon.dirname,'oat.mat'),'oat');

% % Write .nii files
% S                           = [];
% S.oat                       = oat;
% S.stats_fname               = oat.first_level.results_fnames{1};
% S.first_level_contrasts     = [1];
% S.resamp_gridstep           = oat.source_recon.gridstep;
% [statsdir,times,count]      = oat_save_nii_stats(S);

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
