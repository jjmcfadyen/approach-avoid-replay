% Main classification script

clear all
clc

%% Directories

addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults;

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

% Baseline correction on functional localiser epochs
bc = true;
bcwindow = [-.1 0];

% Downsampling
Fs = 100;

% Classification parameters
trainTimes = 0:10:300; % in ms
nT = length(trainTimes);

nulldata = 1; % percent of data to replicate as null
nN = length(nulldata);

nStates = 6;

%% Build classifiers & cross-validate

for s = 1:N
    
    % Load data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_FL_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    data = merged;
    
    % Remove incorrect trials
    behav = readtable(fullfile(dir_behav,subjects{s},[num2str(str2double(subjects{s})) '_fl.csv']));
    idx = strcmp(behav.Acc,'True'); % index of trials to include
    
    data.trialinfo = data.trialinfo(idx,:);
    data.sampleinfo = data.sampleinfo(idx,:);
    data.trial = data.trial(idx);
    data.time = data.time(idx);
    
    % Do baseline correct (if applicable)
    if bc
       
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = bcwindow;
        cfg.continuous = 'no';
    
        data = ft_preprocessing(cfg,data);
        
    end
    
    % plot
    %{
    pdata = cell(1,6);
    for c = 1:6
        cfg = [];
        cfg.trials = find(data.trialinfo==c);
        pdata{c} = ft_timelockanalysis(cfg,data);
    end
    
    figure
    cfg = [];
    cfg.layout = 'CTF275.lay';
    cfg.parameter = 'avg';
    cfg.linecolor = colours(6,'viridis');
    ft_multiplotER(cfg,pdata{:})
    %}
    
    % Create output directory for saving data/labels for the cluster
    dir_save = fullfile(dir_batch,subjects{s});
    if ~exist(dir_save)
        mkdir(dir_save)
    end
    
    % Create classifiers for each time point for this subject
    for t = 1:nT
        for n = 1:nN
            
            tic
            
            disp('===================================')
            disp(['=== SUBJECT ' subjects{s} ', TP = ' num2str(trainTimes(t)) ' ms ==='])
            disp('===================================')
            
            % Organise the data for classification
            [X,Y] = organise_FL_data(data,trainTimes(t),nulldata(n));
            
            % Save for the cluster
            filename = fullfile(dir_save,['data_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(nulldata(n)) '.mat']);
            save(filename,'X','Y');
            
            % Build classifier & cross-validate
            crossvalidate_classifier(filename);

%             % Create job for cluster
%             generate_jobs_classifier(subjects{s},trainTimes(t),nulldata(n));
            
            toc

        end
    end
end

%% Consolidate classifiers

% choose null data
thisnull = 1; % as an index of 'nulldata' variable

% Load and select best lambda (on average)
acc = nan(N,nT,nStates);
for s = 1:N
    for t = 1:nT
            
        % Load 'cv' and 'classifier' variables
        load(fullfile(dir_batch,subjects{s},...
            ['cv_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(nulldata(thisnull)) '.mat']));

        % Average over folds
        cv = squeeze(nanmean(cv));

        % Get best lambda (averaged across conditions)
        [~,best_lambda] = max(nanmean(cv));

        % Store in variable
        acc(s,t,:) = cv(:,best_lambda);

    end
end

% Plot by subject
cmap = colours(N,'plasma');

figure

% plot subjects
for s = 1:N
    d = squeeze(acc(s,:,:))';
    m = mean(d);
    sem = std(d)/size(d,1);
    upper = m+sem;
    lower = m-sem;
    for i = 1:length(sem)
        plot(repmat(trainTimes(i),2,1),[lower(i) upper(i)],'color',cmap(s,:),'linewidth',1); hold on;
    end
    plot(trainTimes,m,'color',cmap(s,:),'linewidth',1.2); hold on; % average over states
end

% plot group average
d = squeeze(mean(acc,3));
gm = mean(d);
gsem = std(d)/sqrt(size(d,1));
gupper = gm+gsem;
glower = gm-gsem;
patch([trainTimes fliplr(trainTimes)],[gupper fliplr(glower)],'k','facealpha',.2,'edgealpha',0); hold on
plot(trainTimes,gm,'k','linewidth',2); hold on

% plot chance level
plot(trainTimes([1 end]),repmat(1/nStates,2,1),'k:'); % hypothetical chance level (i.e., 1/6)

% details
set(gca,'ticklength',[0 0])
xlabel('Training Time (ms)')
ylabel('Classification Accuracy')
title('Classification Accuracy per Subject')

% Plot by state/image
