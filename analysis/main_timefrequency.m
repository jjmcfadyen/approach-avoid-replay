% Time-frequency analysis

clear all
clc

%% Directories

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_meg = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';
dir_scripts = 'D:\2020_RiskyReplay\approach-avoid-replay\';

cd(dir_scripts)

%% Settings

addpath('utils');
addpath('preprocessing')
addpath(genpath('analysis'))

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults

% Replay onsets
lagrange = 20:10:90; % lags at which to look for onsets (in ms)
trainTimes = 0:10:300; % in ms

load('D:\2020_RiskyReplay\data\meg\replay\withoutintercept\optimised_times.mat'); % get best classifier training times per participant

%% Overall power spectra throughout planning period

% Fs = 600;
% 
% PS = cell(N,1);
% T = [];
% for s = 1:N
%     
%     disp('==============================')
%     disp(['=== ' subjects{s} ' ==='])
%     disp('==============================')
%     
%     % Load data
%     load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
%     data = merged;
%     clear merged;
%     
%     nTrls = length(data.trial);
%     
%     % Load behavioural data
%     load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
%     behav = behav.task;
%     
%     % Match with MEG data
%     idx = zeros(nTrls,1);
%     for trl = 1:nTrls
%         if data.trialinfo(trl,1) == 0
%             thispractice = 1;
%             thisblock = 1;
%         else
%             thispractice = 0;
%             thisblock = data.trialinfo(trl,1);
%         end
%         thistrial = data.trialinfo(trl,2);
%         idx(trl,1) = find(behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial);
%     end
%     behav = behav(idx,:);
%     
%     T = [T; behav(:,ismember(behav.Properties.VariableNames,{'Subject','Forced','nV_1','nV_2','P','Choice','RT'}))];
%     
%     % Exclude baseline
%     cfg = [];
%     cfg.toilim = [0 Inf];
%     data = ft_redefinetrial(cfg,data);
%     
%     % Power spectra in planning period of each trial
%     cfg = [];
%     cfg.output = 'pow';
%     cfg.channel = 'MEG';
%     cfg.method = 'mtmfft';
%     cfg.taper = 'boxcar';
%     cfg.foi = 2:150;
%     cfg.keeptrials = 'yes';
%     cfg.pad = 'nextpow2';
%     
%     PS{s} = ft_freqanalysis(cfg,data);
% 
% end
% 
% % Statistics
% % -- forced vs free
% pd = [];
% for s = 1:N
%     
%     thisT = T(contains(T.Subject,subjects{s}),:);
%     
%     % select maximal channels?
%     maxchan = squeeze(mean(mean(PS{s}.powspctrm,3)));
%     maxchan = maxchan > quantile(maxchan,.95);
%     
%     for c = 1:2
%         if c==1
%             idx = find(thisT.Forced~=0);
%         elseif c==2
%             idx = find(thisT.Forced==0);
%         end
%         pd(s,c,:) = squeeze(mean(mean(PS{s}.powspctrm(idx,maxchan,:),2)));
%     end
% end
% 
% x = PS{1}.freq;
% cmap = [0 0 1; 1 0 0];
% figure
% for c = 1:2
%     m = squeeze(mean(pd(:,c,:)))';
%     sem = squeeze(std(pd(:,c,:))/sqrt(size(pd,1)))';
%     upper = m+sem;
%     lower = m-sem;
%     patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
%     plot(x,m,'color',cmap(c,:),'linewidth',1); hold on
% end
% 
% % -- approach vs avoid
% pd = [];
% for s = 1:N
%     
%     thisT = T(contains(T.Subject,subjects{s}),:);
%     
%     % select maximal channels?
%     maxchan = squeeze(mean(mean(PS{s}.powspctrm,3)));
%     maxchan = maxchan > quantile(maxchan,.95);
%     
%     for c = 1:2
%         idx = find(thisT.Choice==c);
%         pd(s,c,:) = squeeze(mean(mean(PS{s}.powspctrm(idx,maxchan,:),2)));
%     end
% end
% 
% x = PS{1}.freq;
% cmap = [255, 170, 39
%     176, 39, 255]/255;
% figure
% for c = 1:2
%     m = squeeze(mean(pd(:,c,:)))';
%     sem = squeeze(std(pd(:,c,:))/sqrt(size(pd,1)))';
%     upper = m+sem;
%     lower = m-sem;
%     patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
%     plot(x,m,'color',cmap(c,:),'linewidth',1); hold on
% end

%% Time-frequency

bigT = [];
for s = 1:N
    
    % Compute here
    directories = struct();
    directories.dir_meg = dir_meg;
    directories.dir_behav = dir_behav;
    directories.dir_classifiers = dir_classifiers;
    directories.dir_save = 'D:\2020_RiskyReplay\data\meg\timefrequency';
    
    thisT = compute_timefrequency(subjects{s},optimised_times(s),directories,true);
    bigT = [thisT; bigT];
    
    % OR make jobs for cluster
    directories = struct();
    directories.dir_meg = '~/Scratch/2020_RiskyReplay/data/meg';
    directories.dir_behav = '~/Scratch/2020_RiskyReplay/data/behav';
    directories.dir_classifiers = '~/Scratch/2020_RiskyReplay/data/meg/classifiers';
    directories.dir_save = '~/Scratch/2020_RiskyReplay/data/meg/timefrequency';
    
    subject = subjects{s};
    optimised_time = optimised_times(s);
    
    filename = [subjects{s} '_tf-data.mat'];
    save(fullfile(dir_batch,filename),'directories','subject','optimised_time');
    generate_jobs_timefrequency(['~/Scratch/2020_RiskyReplay/scripts/',filename]);
    
end
