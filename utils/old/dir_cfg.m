function [sinfo, dinfo] = dir_cfg()

% Set the subjects (including session order & questionnaires)
% Set the directory structure

sinfo = [];
dinfo = [];

%% Subjects

sinfo.subjects = {
    '088674'
    '263098'  % POOR ACCURACY (51.37%)
    '680913'  % POOR ACCURACY (47.22%)
%     '542599'  % ABORTED EXPERIMENT EARLY (BLOCKS 8-10 MISSING) & ISSUES WITH PHOTODIODE DURING FUNCTIONAL LOCALISER
    '383991'
    '707132'
    '396430'
    '521846'
    '015397'
    '663186'
    '097403'  % POOR ACCURACY (60.84%)
    '503608'  % POOR ACCURACY (60.27%)
    '147947'
    '506559'
    };

sinfo.sessions = [0 0 0 1 1 1 2 2 1 2 2 2 2]'; % whether structure A or B was used (0 = older version that had too many catch trials)

sinfo.questionnaires = [ % copied from participantSummary.xlsx on OneDrive
 0.4167 	 0.4625 	 0.6190 	 0.5476 	 0.6667 	 0.5476 	 0.6905 
 0.6000 	 0.6250 	 0.4762 	 0.7619 	 0.2619 	 0.7143 	 0.6429 
 0.7833 	 0.8500 	 0.2619 	 0.3333 	 0.3571 	 0.2143 	 0.4524 
%  0.4667 	 0.7875 	 0.3810 	 0.1429 	 0.2619 	 0.2857 	 0.7857 
 0.3667 	 0.6000 	 0.3333 	 0.1905 	 0.6190 	 0.4286 	 0.5952 
 0.2333 	 0.4000 	 0.2857 	 0.2143 	 0.1429 	 0.3571 	 0.4048 
 0.3833 	 0.7000 	 0.2143 	 0.3810 	 0.4048 	 0.2381 	 0.8571 
 0.6000 	 0.7000 	 0.2143 	 0.3571 	 0.3333 	 0.2143 	 0.7857 
 0.5667 	 0.8375 	 0.2143 	 0.2619 	 0.4286 	 0.2619 	 0.6429 
 0.5833 	 0.6625 	 0.2143 	 0.3095 	 0.5476 	 0.4762 	 0.7857 
 0.7333 	 0.8375 	 0.6190 	 0.6190 	 0.8095 	 0.4762 	 0.8810 
 0.5500 	 0.7875 	 0.1905 	 0.2857 	 0.2619 	 0.4524 	 0.5238 
 0.6167 	 0.7625 	 0.2143 	 0.3095 	 0.1667 	 0.2381 	 0.5952 
 0.4333 	 0.5625 	 0.3810 	 0.4524 	 0.4762 	 0.8571 	 0.7619 
 ]; % intolerance of uncertainty, worry, ethical risk, financial risk, health risk, recreational risk, social risk

maxscores = [60 80 42 42 42 42 42];
sinfo.questionnnaires_orig = round(sinfo.questionnaires .* maxscores);

sinfo.nC = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]; % negator combos

%% Directories

dinfo.home = 'D:\2020_RiskyReplay';

% data
dinfo.data = fullfile(dinfo.home,'data');

dinfo.data_behav = fullfile(dinfo.data,'behav');
dinfo.data_meg = fullfile(dinfo.data,'meg');
dinfo.data_meg_raw = fullfile(dinfo.data_meg,'raw');
dinfo.data_meg_pp1 = fullfile(dinfo.data_meg,'preprocessed_1');
dinfo.data_meg_pp2 = fullfile(dinfo.data_meg,'preprocessed_2');
dinfo.data_meg_pp3 = fullfile(dinfo.data_meg,'preprocessed_3');
dinfo.data_meg_classifiers = fullfile(dinfo.data_meg,'classifiers');

% results
dinfo.results = fullfile(dinfo.home,'results');

dinfo.results_behav = fullfile(dinfo.results,'behav');
dinfo.results_meg = fullfile(dinfo.results,'meg');
dinfo.results_meg_classifiers = fullfile(dinfo.results_meg,'classifiers');

% toolboxes
dinfo.tb = 'D:\Toolboxes';
dinfo.tb_fieldtrip = fullfile(dinfo.tb,'fieldtrip-20191119');
dinfo.tb_spm = fullfile(dinfo.tb,'spm12');
dinfo.tb_osl = fullfile(dinfo.tb,'osl-core');

end