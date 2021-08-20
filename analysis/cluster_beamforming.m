function cluster_preprocess(filename)
% Script to preprocess the data 
% (needs to be run on UNIX based system using OSL toolbox)
% r = run number (double)

%% Directories

dir_osl = '~/Scratch/2020_RiskyReplay/toolboxes/osl-core-master-UNIX/';
addpath(genpath(dir_osl));
osl_startup;

%% Beamforming (all events)

% Load OAT variable
load(filename);

if ~exist(oat.source_recon.dirname)
    mkdir(oat.source_recon.dirname)
end

% Run
oat = osl_check_oat(oat);
oat = osl_run_oat(oat);

end