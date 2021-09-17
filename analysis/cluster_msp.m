function cluster_msp()

%% Set up

dir_spm = '/data/holly-host/jmcfadyen/spm12';
addpath(dir_spm);
spm('defaults','eeg');

subjects = {'012882'
    '013770'
    '015397'
    '018768'
    '027480'
    '088674'
    '097403'
    '147947'
    '220598'
%     '263098'
    '383991'
    '391883'
    '396430'
    '503608'
    '506559'
    '521846'
    '663186'
%     '680913'
    '706130'
    '707132'
    '795776'
    '832746'
    '909566'
    '945092'
    '957849'
    '989833'};
N = length(subjects);

%% Group Inversion

dir_epochs = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/msp';
dir_canonical = fullfile(dir_spm,'canonical');

% Parameters
meshres = 2; % 1 = coarse, 2 = normal, 3 = fine

if meshres==1
    meshunits = 5124;
elseif meshres==2
    meshunits = 8196;
elseif meshres==3
    meshunits = 20484;
end

for i = 1%:2 % 1 = template cortical mesh, 2 = freesurfer colin cortical mesh
    
    % Overwrite template with either the template or a new mesh
    if i==1
        copyfile(fullfile(dir_canonical,'orig',['cortex_' num2str(meshunits) '.surf.gii']),...
            fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
        copyfile(fullfile(dir_canonical,'orig','single_subj_T1.nii'),...
            fullfile(dir_canonical,'single_subj_T1.nii'));
    elseif i==2
        copyfile(fullfile(dir_canonical,'new','canonicalHippAmg.gii'),...
            fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
        copyfile(fullfile(dir_canonical,'new','orig.nii'),...
            fullfile(dir_canonical,'single_subj_T1.nii'));
    end
    
%     % Generate priors
%     if i==2
%         mesh = gifti(fullfile(dir_canonical,['cortex_' num2str(meshunits) '.surf.gii']));
%         meshcoords = spm_eeg_inv_transform_points(mesh.mat,mesh.vertices);
%         % match the vertices with cortical vs hipp vs amg
%     end
    
    filenames = {};
    for s = 1:N
        
        dir_save = fullfile(dir_epochs,subjects{s},['inv' num2str(i)]);
        if ~exist(dir_save)
            mkdir(dir_save)
        end
        
        fname = ['replay-epochs_all_' subjects{s} '.mat'];
        filename = fullfile(dir_save,fname);
        filenames{s,1} = filename;
        
%         D = spm_eeg_load(fullfile(dir_epochs,subjects{s},fname));
%         disp(['Moving ' filename '...'])
%         if ~exist(filename)
%             copy(D,filename);
%         end
        
    end
        
    matlabbatch = {};
    cc = 0;

    % Head model specification
%     cc = cc+1;
%     matlabbatch{cc}.spm.meeg.source.headmodel.D = filenames;
%     matlabbatch{cc}.spm.meeg.source.headmodel.val = 1;
%     matlabbatch{cc}.spm.meeg.source.headmodel.comment = '';
%     matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
%     matlabbatch{cc}.spm.meeg.source.headmodel.meshing.meshres = 2;
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
%     matlabbatch{cc}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
%     matlabbatch{cc}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
%     matlabbatch{cc}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
% 
%     % Source inversion
%     cc = cc+1;
%     matlabbatch{cc}.spm.meeg.source.invert.D = filenames;
%     matlabbatch{cc}.spm.meeg.source.invert.val = 1;
%     matlabbatch{cc}.spm.meeg.source.invert.whatconditions.condlabel = {
%         'aversivereplay_avoid'
%         'rewardingreplay_avoid'
%         'aversivereplay_approach'
%         'rewardingreplay_approach'
%         }';
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.woi = [-100 150];
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
%     matlabbatch{cc}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
%     matlabbatch{cc}.spm.meeg.source.invert.modality = {'MEG'};
% 
%     % Get results
%     cc = cc+1;
%     matlabbatch{cc}.spm.meeg.source.results.D = filenames;
%     matlabbatch{cc}.spm.meeg.source.results.val = 1;
%     matlabbatch{cc}.spm.meeg.source.results.woi = [0 100];
%     matlabbatch{cc}.spm.meeg.source.results.foi = [0 0];
%     matlabbatch{cc}.spm.meeg.source.results.ctype = 'evoked';
%     matlabbatch{cc}.spm.meeg.source.results.space = 1;
%     matlabbatch{cc}.spm.meeg.source.results.format = 'mesh';
%     matlabbatch{cc}.spm.meeg.source.results.smoothing = 8;
%     
%     c = cc+1;
%     matlabbatch{cc}.spm.meeg.source.results.D = filenames;
%     matlabbatch{cc}.spm.meeg.source.results.val = 1;
%     matlabbatch{cc}.spm.meeg.source.results.woi = [0 100];
%     matlabbatch{cc}.spm.meeg.source.results.foi = [0 0];
%     matlabbatch{cc}.spm.meeg.source.results.ctype = 'evoked';
%     matlabbatch{cc}.spm.meeg.source.results.space = 1;
%     matlabbatch{cc}.spm.meeg.source.results.format = 'image';
%     matlabbatch{cc}.spm.meeg.source.results.smoothing = 8;

    % Do 2nd level group analysis
    dir_group = fullfile(dir_epochs,['group_inv' num2str(i)]);
    if ~exist(dir_group)
        mkdir(dir_group)
    end
    if exist(fullfile(dir_group,'SPM.mat'))
        warning(['Deleting ' fullfile(dir_group,'SPM.mat') '...'])
        delete(fullfile(dir_group,'SPM.mat'))
    end
    
    factors = {'choice','pathtype'};
    conditions = { 'aversivereplay_avoid'
        'rewardingreplay_avoid'
        'aversivereplay_approach'
        'rewardingreplay_approach'}; % I think this is the order they appear in but I'm not sure... (this is what it is in D.inv{1}.inverse.trials')
    levels = [2 2;
              2 1
              1 2
              1 1];
    
    cc = cc+1;
    matlabbatch{cc}.spm.stats.factorial_design.dir = {dir_group};
    for f = 1:length(factors)
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).name = factors{f};
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).levels = 2;
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).dept = 0;
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).variance = 1;
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).gmsca = 0;
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.fact(f).ancova = 0;
    end
    for l = 1:size(levels,1)
        matlabbatch{cc}.spm.stats.factorial_design.des.fd.icell(l).levels = levels(l,:);
        for s = 1:N
            matlabbatch{cc}.spm.stats.factorial_design.des.fd.icell(l).scans{s,1} = fullfile(dir_epochs,...
                subjects{s},['inv' num2str(i)],['replay-epochs_all_' subjects{s} '_1_t0_100_f_1.gii']);
        end
    end
    matlabbatch{cc}.spm.stats.factorial_design.des.fd.contrasts = 1;
    matlabbatch{cc}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{cc}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{cc}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{cc}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{cc}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{cc}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{cc}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{cc}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    cc = cc+1;
    matlabbatch{cc}.spm.stats.fmri_est.spmmat = {fullfile(dir_group,'SPM.mat')};
    matlabbatch{cc}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{cc}.spm.stats.fmri_est.method.Classical = 1;
    
    % Run job(s)
    spm_jobman('run',matlabbatch);
end

end