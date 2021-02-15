%% Preprocess data using OSL

% requires raw data .ds folders are placed in the 'data' folder
% requires the s2_data.json and CSV files are placed in the 'data' folder

clear all
close all
clc

%% Set up

% directories
restoredefaultpath;
OSLDIR = '/Users/jmcfadyen/Toolboxes/osl';
addpath(OSLDIR);
cd(fullfile(OSLDIR,'osl-core'));
osl_startup;

dir_home = '/Users/jmcfadyen/Documents/osl_analysis';

%% Organise files

% get list of files to preprocess
dir_data = '/Volumes/SOFTWARE BA/halfprocessed/todo';
addpath(dir_data);
datafiles = dir(fullfile(dir_data,'Bdffp*.mat'));

%% Parameters

params = [];
params.eog = {'UADC001','UADC003','UADC005'};
params.trigger_chan = {'UADC004'};
params.resp_chan = {'UPPT002'};
params.Fs = 100; % frequency to downsample to (in Hz)

params.dir_data = dir_data;

params.epoch = 'decision'; % 'decision', 'transition', 'image', 'outcome'

%% Run

cd(dir_home);

for f = 1:length(datafiles)
    
    fp = strsplit(datafiles(f).name,'_');
    subject = fp{2};
    dtype = fp{3};
    fnum = strsplit(fp{4},'.mat');
    fnum = fnum{1};
    
    if ~exist(fullfile('/Volumes/SOFTWARE BA/processed/',subject,[subject '_' dtype '_' fnum '_ICA.mat']))
        
        % Import
        func_prep(datafiles(f).name,params);

        % ICA
        func_ica(datafiles(f).name,params);
        
    end
    
    % Epoch
    func_epoch(datafiles(f).name,params);
    
end
