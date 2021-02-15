function [endsamples, trigger_table] = preprocessMEG_1(subjects, session, loadEndSamples)

% subjects = cell array of strings - e.g., {'088674','707132'}
% session = 'FL' or 'Task'
% loadEndSamples = true or false (default = false)

% Uses: raw CTF MEG data
% Produces: 
% ---------- MAT files for the point at which to crop the data
% ---------- photodiode trigger onsets & labels

%% Directories

addpath('utils');
[sinfo, dinfo] = dir_cfg();

% Fieldtrip
addpath(dinfo.tb_fieldtrip);
ft_defaults;

%% Parameters

if nargin == 0
    error('No subjects or session type specified');
end
if nargin == 1
    error('No session type specified (FL or Task)');
end
if nargin < 3
    loadEndSamples = false;
end

if ~iscell(subjects)
    subjects = {subjects};
end
        
if ~strcmp(session,'FL') && strcmp(session,'Task')
    error('Wrong session type entered (FL or Task accepted)')
end

if ~islogical(loadEndSamples)
    error('loadEndSamples entered incorrectly - needs to be true or false (default = false)');
end
        
%% Run

for subj = 1:length(subjects)
    
    subject = subjects{subj};
    disp('--------------------------------------------------')
    disp(['--- ' subject '----------------------------------------'])
    disp('--------------------------------------------------')
    disp(['--- (' num2str(round(subj/length(subjects)*100,2)) '%) --------------------------------------'])
    disp('--------------------------------------------------')
    
    [fl,event,in] = pp_cfg(subject,session);
    
    dir_save = fullfile(dinfo.data_meg_pp1,subject);
    if ~exist(dir_save)
        mkdir(dir_save);
    end
    
    %% Manual croppings

    if ~loadEndSamples
        [endSamples,fsave] = pp_findEndSamples(fl,in);
        save( fullfile(dir_save,fsave), 'endSamples');
    else
        fname = dir( fullfile(dinfo.data_meg_pp1,subject,[subject '_' session '_endsamples-*Hz.mat']) );
        load( fullfile(fname.folder,fname.name) );
    end

    %% Identify triggers

    for f = 1:length(fl.data)
        [tmp, fsave] = pp_triggers(fl,event,in,f,endSamples(f));
        if f == 1
            triggerTable = tmp;
        else
            triggerTable = [triggerTable; tmp];
        end
    end
    
    save( fullfile(dinfo.data_meg_pp1,fsave), 'endSamples');
    
end
