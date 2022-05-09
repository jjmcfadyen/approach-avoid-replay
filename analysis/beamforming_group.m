
subjects = [
     {'012882'}
    {'013770'}
    {'015397'}
    {'018768'}
    {'027480'}
    {'088674'}
    {'097403'}
    {'147947'}
    {'220598'}
    {'263098'}
    {'383991'}
    {'391883'}
    {'396430'}
    {'503608'}
    {'506559'}
    {'521846'}
    {'663186'}
    {'680913'}
    {'706130'}
    {'707132'}
    {'795776'}
    {'832746'}
    {'909566'}
    {'945092'}
    {'957849'}
    {'989833'}
    ];

N = length(subjects);

oat = [];

prevoat = load(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/',[subjects{1} '.oat'],...
        'oat.mat'));
oat.source_recon = prevoat.oat.source_recon;

oat.source_recon.dirname = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/group';
if ~exist(oat.source_recon.dirname)
    mkdir(oat.source_recon.dirname);
end

copyfile(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/',[subjects{1} '.oat'],'replay_sub_level_mask.nii.gz'),...
    fullfile([oat.source_recon.dirname '.oat'],'first_level_sub_level_group_level_mask.nii.gz'))

oat.first_level = [];
for s = 1:N
    filelist = dir(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/',[subjects{s} '.oat'],'fsession*_firstlevel.mat'));
    for f = 1:length(filelist)
        oat.first_level.results_fnames{s}(f,1) = {fullfile(filelist(f).folder,filelist(f).name)};
    end
end

oat.subject_level = [];
for s = 1:N
    prevoat = load(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/',[subjects{s} '.oat'],...
        'oat.mat'));
    oat.subject_level.session_index_list{s} = prevoat.oat.subject_level.session_index_list;
    oat.subject_level.results_fnames{s} = fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/beamforming/NOCONT_NEWEPOCHS_anyreplay_4-8Hz/',[subjects{s} '.oat'],...
        'subject1_replay_sub_level.mat');
end
oat.subject_level.subjects_to_do = 1:N;

oat.group_level = [];
oat.group_level.group_design_matrix = ones(1,N);
oat.group_level.group_contrast{1} = [1]; % group average
oat.group_level.subjects_to_do = 1:N;
oat.group_level.time_range = [0 0.1];

oat.to_do = [0 0 0 1];

oat = osl_check_oat(oat);
oat = osl_run_oat(oat);

% view results
S = [];
S.oat = oat;
S.stats_fname = oat.group_level.results_fnames;
S.group_level_contrasts = [1]; % list of first level contrasts to output
S.first_level_contrasts = 1; 
statsdir = oat_save_nii_stats(S);

% cluster permutation
S = [];
S.oat = oat;
S.cluster_stats_thresh = 6;
S.cluster_stats_nperms = 1000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do = [1];
S.group_level_copes_to_do = [1];
S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script = 0;
S.time_range = [0 0.1];
S.time_average = 1;

% Run the permutations
[gstats] = oat_cluster_permutation_testing(S);
