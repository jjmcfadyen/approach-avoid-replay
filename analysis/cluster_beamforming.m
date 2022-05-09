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
cd /data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/batch
load(filename);

% if ~exist(oat.source_recon.dirname)
%     mkdir(oat.source_recon.dirname);
% end

% % Coregister first
% for e = 1:2
% 
%     if e==1
%         fname = 'D_continuous';
%     elseif e==2
%         fname = 'D_epoched';
%     end
% 
%     if isfield(oat.source_recon,fname)
%         for r = 1:length(oat.source_recon.(fname))
%             
%             S = [];
%             S.D = oat.source_recon.(fname){r};
%             
%             if ~isempty(S.D)
%                 D = spm_eeg_load(S.D);
%                 if any(contains(D.chantype,'EEG'))
%                     D = chantype(D,find(contains(D.chantype,'EEG')),'Other');
%                     D.save;
%                 end
%                 
%                 S.mri = []; % empty means to use template
%                 S.useheadshape = 0;
%                 S.use_rhino = 0;
%                 S.forward_meg = 'Single Shell';
%                 S.fid.label.nasion = 'nas';
%                 S.fid.label.lpa = 'lpa';
%                 S.fid.label.rpa = 'rpa';
%                 D = osl_headmodel(S);
%                 D.save;
%             end
%         end
%     end
% end

% Run source reconstruction
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

save(fullfile(oat.source_recon.dirname,'oat.mat'),'oat');

% Run first level
oat.to_do = [0 1 0 0];

nRuns = length(oat.source_recon.D_epoched);
includeRun = ones(1,nRuns);
for r = 1:nRuns
    try
        oat.first_level.sessions_to_do = r; % need to do this because too few trials makes the first level fail (not enough time samples across trials)
        oat = osl_run_oat(oat);
    catch
        includeRun(r) = 0;
    end
end
if all(~includeRun)
    error('No runs included!')
else
    disp(['Runs included: ' num2str(includeRun)]);
end

if nRuns>1
    oat.first_level.sessions_to_do = find(includeRun);
    oat = osl_run_oat(oat);
end

save(fullfile(oat.source_recon.dirname,'oat.mat'),'oat');

% % Run subject level
% oat.to_do = [0 0 1 0];
% oat.subject_level.session_index_list = {find(includeRun)};
% oat = osl_run_oat(oat);
% 
% % Save
% save(fullfile(oat.source_recon.dirname,'oat.mat'),'oat');

% Write .nii files for subject level result
S                           = [];
S.oat                       = oat;
S.stats_fname               = oat.first_level.results_fnames{1}; %[oat.first_level.results_fnames{f} '.mat'];
S.first_level_contrasts     = 1:length(oat.first_level.contrast);
S.resamp_gridstep           = oat.source_recon.gridstep;
[statsdir,times,count]      = oat_save_nii_stats(S);

% Write SPM.mat files
for f = 1:length(oat.first_level.results_fnames)
    S                   = [];
    S.oat               = oat;
    S.stats_fname       = oat.first_level.results_fnames{f};
    [D_tstat, D_cope]   = oat_save_spm_stats(S);
end

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
