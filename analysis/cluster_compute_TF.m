function [avTF,avBTF] = cluster_compute_TF(subject)

addpath('/data/holly-host/jmcfadyen/spm12')
spm('defaults','eeg');

% Load replay epochs
data = spm_eeg_load(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/8_replayepochs-600Hz_nullthresh_bc/merged',['replay-epochs_all_' subject '.mat']));
Fs = data.fsample;
data = ftraw(data);
data.fsample = Fs;

% Make sure there are no NaNs
ridx = nan(length(data.trial),1);
for trl = 1:length(data.trial)
    ridx(trl) = any(isnan(data.trial{trl}(:)));
end
disp(['Removing ' num2str(sum(ridx)) ' trials for missing data in epoch'])

if any(ridx)
    cfg = [];
    cfg.trials = find(~ridx);
    data = ft_selectdata(cfg,data);
end

% Artifically pad the data (function in Fieldtrip isn't working for some reason)
padwindow = [-1 1]; % around existing window
padded = data;
for trl = 1:length(data.trial)
    padded.trial{trl} = ft_preproc_padding(data.trial{trl}, 'zero', abs(round(padwindow(1)*data.fsample)), abs(round(padwindow(end)*data.fsample)));
    padded.time{trl} = linspace(data.time{trl}(1)+padwindow(1),data.time{trl}(end)+padwindow(end),size(padded.trial{trl},2));
end

% Wavelet analysis
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp([num2str(length(padded.trial)) ' TRIALS DETECTED'])
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cfg = [];
data = ft_selectdata(cfg,padded);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'wavelet';
cfg.width        = 4;
cfg.foi          = 1:150;
cfg.toi          = -0.1:.01:.3;

TF = ft_freqanalysis(cfg, data);

% Save
dir_save = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg/timefrequency';
save(fullfile(dir_save,['tf_' subject '.mat']),'TF');

end