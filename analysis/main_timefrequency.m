% Time-frequency analysis

clear all
clc

%% Directories

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_meg = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';
dir_scripts = 'D:\2020_RiskyReplay\approach-avoid-replay\';
dir_classifiers = 'D:\2020_RiskyReplay\data\meg\classifiers';

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

wavelets = [3 5 7 9 11];

for s = 1:N
    
%     % Compute here
%     directories = struct();
%     directories.dir_meg = dir_meg;
%     directories.dir_behav = dir_behav;
%     directories.dir_save = 'D:\2020_RiskyReplay\data\meg\timefrequency';
%     
%     for w = wavelets
%         compute_timefrequency(subjects{s},directories,w);
%     end
    
    % make jobs for cluster
    directories = struct();
    directories.dir_meg = '~/Scratch/2020_RiskyReplay/data/meg';
    directories.dir_behav = '~/Scratch/2020_RiskyReplay/data/behav';
    directories.dir_save = '~/Scratch/2020_RiskyReplay/data/meg/timefrequency';
    
    subject = subjects{s};
    for w = wavelets
        waveletwidth = w;
        filename = [subjects{s} '_w' num2str(w) '_tf-data.mat'];
        save(fullfile(dir_batch,filename),'directories','subject','waveletwidth');
        generate_jobs_timefrequency(['~/Scratch/2020_RiskyReplay/scripts/',filename],'holly');
    end
end

% writetable(bigT,'D:\2020_RiskyReplay\data\meg\timefrequency\tf_table.csv')
% 
% % LME on big table
% bigT.Choice = bigT.Choice+10;
% bigT.Choice(bigT.Choice==11) = 1; % approach
% bigT.Choice(bigT.Choice==12) = 0; % avoid
% bigT.Subject = categorical(bigT.Subject,unique(bigT.Subject),unique(bigT.Subject));
% 
% lme = fitglme(bigT,'Choice~Theta+Alpha+Beta+LowGamma+HighGamma+RT+(1|Subject)','distribution','binomial')

%% Visualise the replay onsets

avreplay = [];
for s = 1:N
    
    disp('==============================')
    disp(['=== ' subjects{s} ' ==================='])
    disp('==============================')

    % Load merged data file
    load(fullfile(directories.dir_meg,['7_merged_ds-600Hz'],[subjects{s} '_task_600Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;
    
    nTrls = length(data.trial);
    
    % Load behavioural data
    load(fullfile(directories.dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
    behav = behav.task;
    
    % Match with MEG data
    idx = zeros(nTrls,1);
    for trl = 1:nTrls
        if data.trialinfo(trl,1) == 0
            thispractice = 1;
            thisblock = 1;
        else
            thispractice = 0;
            thisblock = data.trialinfo(trl,1);
        end
        thistrial = data.trialinfo(trl,2);
        idx(trl,1) = find(behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial);
    end
    behav = behav(idx,:);
    
    % Load time-frequency
    load(['D:\2020_RiskyReplay\data\meg\timefrequency\' subjects{s} '_tf.mat']);
    
    % Average across choice types
%     figure
    for c = 1:2
        
        % merge
        idx = find(behav.Choice==c & behav.Forced==0);
        tmp = [];
        for i = 1:length(idx)
            cc = size(tmp,1);
            if ~any(isnan(squash(rbTF{idx(i)})))
                try
                    tmp((cc+1):(cc+size(rbTF{idx(i)},1)),:,:) = rbTF{idx(i)};
                catch
                end
            end
        end
%         x = linspace(-150,150,size(tmp,3));
        avreplay(s,c,:,:) = squeeze(mean(tmp));
%         
%         subplot(1,2,c)
%         imagesc(squeeze(mean(tmp)),'XData',x,'YData',freq{1});
%         colormap(colours(100,'inferno'))
%         view(180,90);
%         set(gca,'xdir','reverse')
%         hold on
%         plot([0 0],freq{1}([1 end]),'w:')
%         if c==1
%             title('Approach')
%         elseif c==2
%             title('Avoid')
%             colorbar
%             tmp = squash(squeeze(avreplay(s,:,:,:)));
%             thisclim = [-max(abs(tmp)) max(abs(tmp))];
%             subplot(1,2,1)
%             caxis(thisclim)
%             subplot(1,2,2)
%             caxis(thisclim)
%         end
    end
%     drawnow
%     sgtitle(subjects{s})
    
end

% plot average
x = linspace(-100,150,size(avreplay,4));

figure
for c = 1:2
    subplot(1,2,c)
    imagesc(squeeze(mean(avreplay(:,c,:,:))),'XData',x,'YData',freq{1});
    colormap(colours(100,'inferno'))
    view(180,90);
    set(gca,'xdir','reverse')
    hold on
    plot([0 0],freq{1}([1 end]),'w:')
    if c==1
        title('Approach')
    elseif c==2
        title('Avoid')
        tmp = squash(squeeze(avreplay(s,:,:,:)));
        thisclim = [-max(abs(tmp)) max(abs(tmp))];
        subplot(1,2,1)
        caxis(thisclim)
        subplot(1,2,2)
        caxis(thisclim)
    end
end
