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

% Replay onsets
lagrange = 20:10:90; % lags at which to look for onsets (in ms)
trainTimes = 0:10:300; % in ms

load('D:\2020_RiskyReplay\data\meg\replay\withoutintercept\optimised_times.mat'); % get best classifier training times per participant

%% See if replay onsets are anticorrelated across trial time

T = [];
for s = 1:N
    
    subject = subjects{s};
    
    disp('==============================')
    disp(['=== ' subject ' ==='])
    disp('==============================')

    % Load data
    load(fullfile(dir_meg,['7_merged_ds-100Hz'],[subject '_task_100Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;

    nTrls = length(data.trial);

    % Load behavioural data
    load(fullfile(dir_behav,subject,[subject '_parsedBehav.mat']))
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

    % Get classifier info for this subject
    [~,lambdas] = get_bestLambdas(subject,trainTimes,1);
    load(fullfile(dir_classifiers,subject,['classifier_' subject '_t' num2str(optimised_times(s)) '_n1.mat'])); % loads 'classifier'
    classifier.betas = squeeze(classifier.betas(:,:,lambdas(trainTimes==optimised_times(s))));
    classifier.intercepts = squeeze(classifier.intercepts(:,lambdas(trainTimes==optimised_times(s))));

    trialdur = nan(nTrls,1);
    for trl = 1:nTrls
        trialdur(trl,1) = length(data.time{trl});
    end
    
    pvals = [];
    rhos = [];
    meanamp = nan(nTrls,2,max(trialdur));
    pathfreq = [];
    pathphase = [];
    for trl = 1:nTrls

        disp(['Trial ' num2str(trl) ' of ' num2str(nTrls) '...'])
        
        cfg = [];
        cfg.trials = trl;
        thistrial = ft_selectdata(cfg,data); % 100 hz data for replay

        [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
        
        seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
        path1 = mean(seqevidence(:,1:2),2);
        path2 = mean(seqevidence(:,3:4),2);
        
        [r,p] = corr(path1,path2);
        rhos(trl,1) = r;
        pvals(trl,1) = p;
        
        if behav.Forced(trl)==0 && ((behav.nV_1(trl)>0 && behav.nV_2(trl)<0) || (behav.nV_1(trl)<0 && behav.nV_2(trl)>0))
            if behav.nV_1(trl) > behav.nV_2(trl)
                meanamp(trl,1,1:length(seqevidence)) = path1;
                meanamp(trl,2,1:length(seqevidence)) = path2;
            else
                meanamp(trl,1,1:length(seqevidence)) = path2;
                meanamp(trl,2,1:length(seqevidence)) = path1;
            end
        end
        
        [freq,amp,phase] = getPhase(thistrial.time{1},path1);
        pathfreq(trl,1) = round(freq(amp==max(amp)));
        pathphase(trl,1) = phase(amp==max(amp));
        
        [freq,amp,phase] = getPhase(thistrial.time{1},path2);
        pathfreq(trl,2) = round(freq(amp==max(amp)));
        pathphase(trl,2) = phase(amp==max(amp));
        
    end
    
    % Show correlation between paths
    figure
    subplot(2,1,1)
    histogram(rhos); title('Rho')
    subplot(2,1,2)
    histogram(pvals); title('P values')
    sgtitle(subjects{s})
    
    % Show timecourse averaged across trials (display time up to shortest trial)
    cmap = [0, 223, 115
        245, 0, 82]/255;
    
    figure;
    x = data.time{find(trialdur==min(trialdur),1,'first')};
    for p = 1:2
        plot(x,squeeze(nanmean(meanamp(:,p,1:min(trialdur))))','color',cmap(p,:),'linewidth',1.1); hold on
    end
    title(['~' num2str(round(mean(pathfreq(:)))) ' Hz'])
    
    legend({[num2str(round(mean(pathfreq(:,1)))) ' Hz, ' num2str(round(mean(pathphase(:,1)),2))],...
        [num2str(round(mean(pathfreq(:,2)))) ' Hz, ' num2str(round(mean(pathphase(:,2)),2))]})
    
    % Make table
    thisT = behav(:,ismember(behav.Properties.VariableNames,...
        {'Practice','Block','Trial','Forced','ExpTrial','P','nV_1','nV_2','EV','Choice','Acc','RT','Subject'}));
    thisT.Replay_correlation = rhos;
    thisT.Replay_corrpvals = pvals;
    thisT.Rewarding_freq = pathfreq(:,1);
    thisT.Aversive_freq = pathfreq(:,2);
    thisT.Rewarding_phase = pathphase(:,1);
    thisT.Aversive_phase = pathphase(:,2);
    
    T = [T; thisT];
    
end




