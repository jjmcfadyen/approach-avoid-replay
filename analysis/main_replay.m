% Main replay script

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

% Downsampling
Fs = 100;

% Classification parameters
trainTimes = 0:10:300; % in ms
nT = length(trainTimes);
thisnull = 1;
nStates = 6;

% Optimisation parameters
top_percentile = .2; % after sorting the classifiers by accuracy, within which top percentile to optimise for sequenceness?
group_best_time = 120; % from group average cross-validation accuracy

% Sequenceness parameters
U = generate_nullperms([]); 

%% Select classifier training time by optimising for general sequenceness

optimal_trainingtimes = nan(N,1);
allreplay = [];
for s = 1:N
   
    %% Determine best classifier training times
    CV = [];
    lambdas = [];
    for t = 1:nT
        
        % Load cross-validation accuracy
        load(fullfile(dir_batch,subjects{s},...
            ['cv_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(thisnull) '.mat']));
        
        % Average over folds
        cv = squeeze(nanmean(cv));

        % Get best lambda (averaged across conditions)
        [~,best_lambda] = max(nanmean(cv));
        lambdas(t) = best_lambda;

        % Average over states
        CV(t) = mean(cv(:,best_lambda)); 
        
    end
    
    % Find best training time
    best_tp = find(CV > 1.5/6); % 1.5 x chance level
    best_tp = unique([find(trainTimes==group_best_time) best_tp]); % add group maximum accuracy, if not already included
    
    %% Get replay for each of the best classifiers
    
    % Restrict the null permutations to those with unique transitions in any order
    % (we're going to average across paths 1 and 2, so the order doesn't matter)
    opts = [];
    opts.ignoreOrder = true;
    uIdx = restrict_nullperms(U,opts);
    thisu = U(uIdx,:);

    nPerm = size(thisu,1);
    
    % Get task data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged;
    
    if ~isfield(data,'fsample')
        data.fsample = Fs;
    end
    
    % Get replay for each classifier training time
    dir_classifier = fullfile(dir_batch,subjects{s}); % this is where the classifier data is stored
    
    replay = nan(length(best_tp),length(data.trial),size(thisu,1),3,60);
    for t = 1:length(best_tp)
       
        % build classifier
        filename = fullfile(dir_classifier,['data_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(thisnull) '.mat']);
        classifier = build_classifier(filename);
        classifier.betas = squeeze(classifier.betas(:,:,lambdas(best_tp(t))));
        classifier.intercepts = classifier.intercepts(:,lambdas(best_tp(t)));
        
        % get replay
        thisreplay = compute_replay(data,classifier,thisu);
        replay(t,:,:,:,:) = squeeze(mean(thisreplay,4));
    end
    
    legend(strsplit(num2str(best_tp)))
    
    % get difference between replay and null threshold
    avreplay = squeeze(mean(replay(:,:,:,3,:),2));
    nulldiff = [];
    for t = 1:length(best_tp)
        np = squeeze(avreplay(t,2:end,:));
        abs_np = unique(round(abs(np),4),'rows');
        npThresh = quantile(max(abs_np,[],2),.975);
        nulldiff(t,:) = abs(squeeze(avreplay(t,1,:))) - npThresh;
    end
    
    % get maximum peaks per training time
    [peaks,peaklocs] = max(nulldiff,[],2);
    idx = find(peaklocs > (median(peaklocs)-round(std(peaklocs)/2)) & peaklocs < (median(peaklocs)+round(std(peaklocs)/2))); % look in window around median peak location 
    
    peaks = peaks(idx,:);
    
    startdist = peaks - nulldiff(idx,1);
    criteria = [startdist peaks];
    
    best_of_both = [max(startdist)/2 max(peaks)]; % we want something with the biggest difference between start and peak (indicating low autocorrelation) and the maximum peak
    
    dist = sqrt(sum(bsxfun(@minus, criteria, best_of_both).^2,2));
    
    closest = idx(dist==min(dist));
    thisbestcv = find(best_tp == find(CV==max(CV)));
    overallbestcv = find(best_tp == find(trainTimes==group_best_time));
    
    % plot
    figure
    set(gcf,'position',[7 558 1911 420])
    cmap = colours(4,'plasma');
    for i = 1:3
        if i==1
            thisidx = overallbestcv;
        elseif i==2
            thisidx = thisbestcv;
        elseif i==3
            thisidx = closest;
        end
        for g = 1:3
            subplot(i,3,g)
            opts = [];
            opts.g = g;
            opts.avCol = cmap(i,:);
            if g < 3
                opts.subtractNull = true;
            else
                opts.subtractNull = false;
            end
            d = squeeze(replay(closest,:,:,g,:));
            ss_plot(d,opts); hold on
            drawnow
        end
        if i==1
            sgtitle(['Group: ' num2str(group_best_time) 'ms'])
        elseif i==2
            sgtitle(['Individual: ' num2str(trainTimes(best_tp(thisbestcv))) 'ms'])
        elseif i==3
            sgtitle(['Group: ' num2str(trainTimes(best_tp(closest))) 'ms'])
        end
    end

    % Store original and optimised (averaging over trials)
    allreplay(s,1,:,:,:) = squeeze(nanmean(replay(overallbestcv,:,:,:,:),2));
    allreplay(s,2,:,:,:) = squeeze(nanmean(replay(thisbestcv,:,:,:,:),2));
    allreplay(s,3,:,:,:) = squeeze(nanmean(replay(closest,:,:,:,:),2));
    
    optimal_trainingtimes(s,:) = best_tp(closest);
    
    clear replay
    clear nulldiff
    
end

% Plot the group average with the winning classifer vs replay-optimised classifier

cmap = [0 0 1; 1 0 0]; % blue for original, red for optimised

figure
opts = [];
for i = 1:2
    for g = 1:3
        
        d = squeeze(allreplay(:,i,:,g,:));
        
        subplot(1,3,g)
        opts.g = g;
        opts.avCol = cmap(i,:);
        ss_plot(d,opts);
    end
end
    
