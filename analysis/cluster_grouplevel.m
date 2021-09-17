function cluster_grouplevel(subject)
% Script to preprocess the data
% (needs to be run on UNIX based system using OSL toolbox)
% r = run number (double)

%% Directories

dir_osl = '/data/holly-host/jmcfadyen/osl-core-master-UNIX';
addpath(genpath(dir_osl));
osl_startup;

subjects = {'012882'
    '013770'
    '015397'
    '018768'
    '027480'
    '088674'
    '097403'
    '147947'
    '220598'
    '263098'
    '383991'
    '391883'
    '396430'
    '503608'
    '506559'
    '521846'
    '663186'
    '680913'
    '706130'
    '707132'
    '795776'
    '832746'
    '909566'
    '945092'
    '957849'
    '989833'};
N = length(subjects);

%% Group-level (beamforming)

oat = [];

% Subjects
for s = 1:N

    r = -1;
    while true
        r = r+1;
        thisfile = fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming',[subjects{s} '.oat'],...
                    'subject1_wholebrain_first_level_sub_level.mat');
        if exist(thisfile)
            break
        end
        if r>11
            warning(['No file found for ' subjects{s}]);
            break
        end
    end

    if s==1
        disp(['Loading OAT for ' subjects{s} ' to use as template...'])
        oatin = load(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming',[subjects{s} '.oat'],...
                    'oat.mat'));
        oat = oatin.oat;
    end
    
%     oat.subject_level.results_fnames{s} = thisfile;
%     oat.subject_level.session_index_list{s} = oat.subject_level.session_index_list{1};
%     oat.subject_level.subjects_to_do(s) = s;

end

oat.group_level.name='group_level';
oat.group_level.subjects_to_do = 1:N;

% Spatial and temporal averaging options
oat.group_level.time_range          = [-0.1 0.1];
oat.group_level.space_average       = 0;
oat.group_level.time_average        = 0;
oat.group_level.time_smooth_std     = 0; % secs
oat.group_level.use_tstat           = 0;

% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm                 = 0; % mm
oat.group_level.group_varcope_time_smooth_std       = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm   = 100; % smooths the variance of the group copes. It is recommended to do this.

% Set up design matrix and contrasts
oat.group_level.group_design_matrix     = ones(1,length(oat.group_level.subjects_to_do));
oat.group_level.group_contrast          = [];
oat.group_level.group_contrast{1}       = [1]';
oat.group_level.group_contrast_name     = {};
oat.group_level.group_contrast_name{1}  = 'mean';

% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do         = [1]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do       = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do       = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes       = 0;
oat.group_level.report.show_lower_level_cope_maps   = 0;

oat.to_do = [0 0 0 1];

oat = osl_check_oat(oat);
oat = osl_run_oat(oat);

% Plots
current_dir = fileparts(oat.fname);

S = [];
S.oat = oat;
S.stats_fname = oat.group_level.results_fnames;
S.subject_level_contrasts = [1];
S.first_level_contrasts = [1];
S.group_level_contrasts = [1];
S.resamp_gridstep = 5 ;
[statsdir,times] = oat_save_nii_stats(S);

% ROIs
dir_mask = '/data/holly-host/jmcfadyen/2020_RiskyReplay/mri';
masks = {'leftHipp.nii','rightHipp.nii','leftAmg.nii','rightAmg.nii'};
for m = 1:length(masks)
    S=[];
    S.oat=oat;
    S.stats_fname=oat.group_level.results_fnames;
    S.mask_fname=fullfile(dir_mask,masks{m});
    [stats,times,mni_coords_used]=oat_output_roi_stats(S);

    S=[];
    S.stats=stats;
    S.oat=oat;
    S.first_level_cons_to_do = [1];
    S.group_level_cons_to_do = [1];
    [H] = oat_plot_vox_stats(S);

    saveas(H,fullfile(current_dir,'plots','wholebrain_first_level_sub_level_group_level',....
        [masks{m}(1:end-4) '_ROI.png']));
end

% % MNI coordinates
% mni_coords = [
%     18 -12 -27 % from Yunzhe's 2019 Cell paper (hippocampus)
%     16 -11 -21 % from Yunzhe's 2019 Cell paper (hippocampus)
%     20 -97 -13 % from Yunzhe's 2019 Cell paper (visual cortex)
% ];
% 
% for m = 1:size(mni_coords,1)
%     S2=[];
%     S2.oat=oat;
%     S2.stats=oat.group_level.results_fnames;
%     S2.vox_coord=mni_coords(m,:);
%     S2.first_level_cons_to_do = [1];
%     S2.group_level_cons_to_do = [1];
%     H = oat_plot_vox_stats(S2);
%     saveas(H,fullfile(current_dir,'plots','wholebrain_first_level_sub_level_group_level',...
%           ['MNI_' num2str(mni_coords(m,1)) '_' num2str(mni_coords(m,2)) '_' num2str(mni_coords(m,3)) '.png']));
% end

% Permutation testings
S                                       = [];
S.oat                                   = oat;
S.cluster_stats_thresh                  = 6;
S.cluster_stats_nperms                  = 1000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do               = [1];
S.group_level_copes_to_do               = [1];
S.group_varcope_spatial_smooth_fwhm     = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script                  = 0;
S.time_range                            = [0 0.1];
S.time_average                          = 1;

% Run the permutations (uses FSL's 'randomise' function)
[gstats] = oat_cluster_permutation_testing(S);

% % View permutation stats
% con = S.first_level_copes_to_do(1);
% tstat = fullfile(gstats.dir,['tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz']);
% clus_tstat = fullfile(gstats.dir,['clustere_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz']);
% corr_clus_tstat = fullfile(gstats.dir,['clustere_corrp_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz']);

% Move subject maps to do group-level statistics in SPM on work PC
dir_save = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/subject_images';
if ~exist(dir_save)
    mkdir(dir_save)
end

for s = 1:N
   
    S                       = [];
    S.oat                   = load(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/',...
                                            [subjects{s} '.oat'],'subject1_wholebrain_first_level_sub_level.mat'));
    S.oat                   = S.oat.oat;
    S.stats_fname           = S.oat.first_level.results_fnames{1};
    S.first_level_contrasts = [1];
    S.resamp_gridstep       = S.oat.source_recon.gridstep;
    S.time_range            = [0 0.1];
    [statsdir,times,count]  = oat_save_nii_stats(S);
    
    
    
end

end
