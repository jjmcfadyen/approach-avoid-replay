function func_ica(datafile,params)

tic

[subject,session,event,in] = pp_osl_cfg(datafile);
disp(['----- FILENAME: ' datafile])
disp(['----- SUBJECT:  ' subject])
disp(['----- SESSION:  ' session])

% update output directory
params.dir_output = fullfile('/Volumes/SOFTWARE BA 1/processed',subject);
addpath(params.dir_output);

%% ICA

rng(str2double(subject));

opt3 = [];

% files
opt3.spm_files = {fullfile(params.dir_data,datafile)};
opt3.dirname = fullfile(params.dir_output,'opt3');
opt3.datatype = 'ctf';

% ica
opt3.africa.do = 1;
opt3.africa.ident.max_num_artefact_comps = 20;
opt3.africa.ident.artefact_chans = {'EOG'};
opt3.africa.precompute_topos = false;
opt.africa.ident.do_kurt = true;
opt3.africa.ident.mains_kurt_thresh = 0.5;
opt3.africa.ident.do_cardiac = true;
opt3.africa.ident.do_plots = false;
opt3.africa.ident.do_mains = false;

% switch off
opt3.downsample.do = 0;
opt3.highpass.do = 0;
opt3.bad_segments.do = 0;
opt3.epoch.do = 0;
opt3.outliers.do = 0;
opt3.coreg.do = 0;

% run
opt3 = osl_run_opt(opt3);

%% Save & clean up

fparts = strsplit(datafile,'_');
dtype = fparts{3}; % 'FL' or 'Task'
fnum = strsplit(fparts{end},'.mat');
fnum = str2num(fnum{1});

copyfile( fullfile(params.dir_data,[opt3.results.spm_files_basenames{1} '.mat']),...
          fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.mat']));
delete(fullfile(params.dir_data,[opt3.results.spm_files_basenames{1} '.mat']));
      
copyfile( fullfile(params.dir_data,[opt3.results.spm_files_basenames{1} '.dat']),...
          fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.dat']));
delete(fullfile(params.dir_data,[opt3.results.spm_files_basenames{1} '.dat']));

% update filepath of datafile
load(fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.mat']));
D.path = params.dir_output;
D.fname = [subject '_' dtype '_' num2str(fnum) '_ICA.mat'];
D.data.fname = fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.dat']);

% save
save(fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.mat']),'D');

end