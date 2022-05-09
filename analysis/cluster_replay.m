function cluster_replay(subject,timepoint,nulldata,lambda)

Fs = 600;

dir_classifiers = '/data/holly-host/jmcfadyen/2020_RiskyReplay/data/meg/classifiers/';
dir_replay = fullfile(['/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/replay/withoutintercept_' num2str(Fs) 'Hz/'],subject);
dir_data = ['/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/7_merged_ds-' num2str(Fs) 'Hz/'];

addpath('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/utils')

if ~exist(dir_replay)
    mkdir(dir_replay)
end

tmp = load(fullfile(dir_data,[subject '_task_' num2str(Fs) 'Hz.mat']));
data = tmp.merged;
clear tmp;

if ~isfield(data,'fsample')
    data.fsample = Fs;
end

U = generate_nullperms([]); 

% % build & save classifier
% filename = fullfile(dir_classifiers,subject,['data_' subject '_t' num2str(timepoint) '_n' num2str(nulldata) '.mat']);
% classifier = build_classifier(filename);
% save(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(timepoint) '_n' num2str(nulldata) '.mat']),'classifier');

load(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(timepoint) '_n' num2str(nulldata) '.mat']));

% get replay using best lambdas from classifier & save
classifier.betas = squeeze(classifier.betas(:,:,lambda));
classifier.intercepts = classifier.intercepts(:,lambda);
replay = compute_replay(data,classifier,U);

save(fullfile(dir_replay,['replay_' subject '_t' num2str(timepoint) '_n' num2str(nulldata) '.mat']),'replay');

end
