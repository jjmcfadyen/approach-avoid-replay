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
dir_classifiers = 'D:\2020_RiskyReplay\data\meg\classifiers';

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

makePlots = false;
onCluster = true;

if makePlots
    allstatenames = {'baby','bicycle','bowtie','backpack','car','cat','cupcake','house','zebra','toothbrush','hourglass','lamp'};
    sdata = cell(N,length(allstatenames));
end

tc = nan(N,6);
stateorder = cell(N,6); % the order of images in 'allstatenames' for this subject
for s = 2:N
    
    % Load data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_FL_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged
    
    % Remove incorrect trials
    behav = readtable(fullfile(dir_behav,subjects{s},[num2str(str2double(subjects{s})) '_fl.csv']));
    idx = strcmp(behav.Acc,'True'); % index of trials to include
    
    % Remove outlier RT trials
    idx = idx & abs(zscore(behav.RT)) < 5;
    
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
    
    % Log states
    cc = 0;
    for p = 1:2
        for st = 1:3
            thisidx = idx & behav.Path==(p-1) & behav.State==(st-1);
            cc = cc+1;
            tc(s,cc) = sum(thisidx);
            stateorder{s,cc} = behav.Image(find(idx,1,'first'));
        end
    end
    
    % Make plots
    if makePlots
        for st = 1:length(allstatenames)
            thisidx = ismember(behav.Image(idx),allstatenames{st});
            if any(thisidx)
                cfg = [];
                cfg.trials = find(thisidx);
                sdata{s,st} = ft_timelockanalysis(cfg,data);
            end
        end
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
            
            disp('===================================')
            disp(['=== SUBJECT ' subjects{s} ', TP = ' num2str(trainTimes(t)) ' ms ==='])
            disp('===================================')
            
            % Organise the data for classification
            [X,Y] = organise_FL_data(data,trainTimes(t),nulldata(n));
            
            % Save for the cluster
            filename = fullfile(dir_save,['data_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(nulldata(n)) '.mat']);
            save(filename,'X','Y');

            if ~onCluster

                % Build classifier & cross-validate
                classifier = build_classifier(filename);
                save(fullfile('D:\2020_RiskyReplay\data\meg\classifiers',subjects{s},...
                    ['classifier_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(nulldata(n)) '.mat']),'classifier');
                crossvalidate_classifier(filename);
            else

                % Create job for cluster
                generate_jobs_classifier(subjects{s},trainTimes(t),nulldata(n),'holly');
            end
        end
    end
end

% plot group average 
if makePlots
    
    adata = cell(1,12);
    for st = 1:12
        adata{st} = ft_timelockgrandaverage([],sdata{~cellfun('isempty',sdata(:,st)),st});
    end
    gadata = ft_timelockgrandaverage([],adata{:});
    
%     figure
%     tp = [0:0.05:0.3];
%     for t = 1:length(tp)
%         subplot(1,length(tp),t)
%         cfg = [];
%         cfg.layout = 'CTF275.lay';
%         cfg.marker = 'off';
%         cfg.style = 'straight';
%         cfg.comment = 'no';
%         cfg.xlim = repmat(tp(t),1,2);
%         cfg.zlim = [min(gadata.avg(:)) max(gadata.avg(:))];
%         cfg.colormap = colours(100,'viridis');
%         ft_topoplotER(cfg,gadata);
%         title([num2str(tp(t)*1000) ' ms'])
%     end

    thischan = ismember(adata{1}.label,'MRO33');

    cmap = colours(12,'viridis');
    maxval = nan(12,1);
    for st = 1:12
        maxval(st,1) = max(adata{st}.avg(thischan,:));
    end
    [~,sortorder] = sort(maxval);
    
    figure
    for st = 1:12
        x = adata{sortorder(st)}.time;
        m = adata{sortorder(st)}.avg(thischan,:);
        sem = adata{sortorder(st)}.sem(thischan,:);
        upper = m+sem;
        lower = m-sem;
        patch([x fliplr(x)],[upper fliplr(lower)],cmap(st,:),'facealpha',.2,'edgealpha',0,'handlevisibility','off'); hold on
        plot(adata{st}.time,m,'color',cmap(st,:),'linewidth',1.2); hold on
    end
    L = legend(allstatenames(sortorder));
    L.Title.String = 'States';
    set(gca,'ticklength',[0 0])
    xlabel('Time (seconds')
    ylabel('Amplitude')

end

%% Consolidate classifiers

% choose null data
thisnull = 1; % as an index of 'nulldata' variable

% Load and select best lambda (on average)
acc = nan(N,nT,nStates);
best_lambdas = nan(N,nT);
for s = 1:N
    for t = 1:nT
            
        % Load 'cv' and 'classifier' variables
        load(fullfile(dir_classifiers,subjects{s},...
            ['cv_' subjects{s} '_t' num2str(trainTimes(t)) '_n' num2str(nulldata(thisnull)) '.mat']));

        % Average over folds
        cv = squeeze(nanmean(cv));

        % Get best lambda (averaged across conditions)
        [~,best_lambda] = max(nanmean(cv));
        best_lambdas(s,t) = best_lambda;

        % Store in variable
        acc(s,t,:) = cv(:,best_lambda);

    end
end

%% Plot by subject

cmap = colours(N,'plasma');
[sorted,sortidx] = sort(max(mean(acc,3),[],2));

figure

% plot subjects
for s = 1:N
    d = squeeze(acc(sortidx(s),:,:))';
    m = mean(d);
%     sem = std(d)/size(d,1);
%     upper = m+sem;
%     lower = m-sem;
    upper = max(acc(sortidx(s),:,:),[],3);
    lower = min(acc(sortidx(s),:,:),[],3);
    for i = 1:nT
        plot(repmat(trainTimes(i),2,1),[lower(i) upper(i)],'color',cmap(s,:),'linewidth',1,'handlevisibility','off'); hold on;
    end
    plot(trainTimes,m,'color',cmap(s,:),'linewidth',1.2); hold on; % average over states
end

% plot group average
d = squeeze(mean(acc,3));
gm = mean(d);
gsem = std(d)/sqrt(size(d,1));
gupper = gm+gsem;
glower = gm-gsem;
patch([trainTimes fliplr(trainTimes)],[gupper fliplr(glower)],'k','facealpha',.2,'edgealpha',0,'handlevisibility','off'); hold on
plot(trainTimes,gm,'k','linewidth',2); hold on

% plot chance level
plot(trainTimes([1 end]),repmat(1/nStates,2,1),'k:','handlevisibility','off'); % hypothetical chance level (i.e., 1/6)

% details
set(gca,'ticklength',[0 0])
xlabel('Training Time (ms)')
ylabel('Classification Accuracy')
title('Classification Accuracy per Subject')
colorbar;

%% Plot best subject's beta maps

besttime = 120; %trainTimes(mean(acc(s,:,:),3)==max(mean(acc(s,:,:),3)));

[sorted,sortidx] = sort(max(mean(acc,3),[],2));
s = sortidx(end-4);

layout = load('D:\Toolboxes\fieldtrip-20191119\template\layout\CTF275_helmet.mat');
pos = layout.lay.pos;

for s = 1:N
    
    % Load classifier
    load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(besttime) '_n1.mat']));

    % Load data and replace with betas
    load(fullfile('D:\2020_RiskyReplay\data\meg\7_merged_ds-100Hz',[subjects{s} '_FL_100Hz.mat']));
    data = merged;
    avg = ft_timelockanalysis([],data);

    % Make sensor plots of beta weights per state
    thesebetas = [];
    for st = 1:6
        thesebetas(st,:) = squeeze(classifier.betas(st,:,best_lambdas(s,trainTimes==besttime)));
    end
    
    cmap = colours(100,'redblue');
    chanidx = find(ismember(layout.lay.label,avg.label));

    % beta sensor maps
    figure
    set(gcf,'position',[1 -149 1920 973])
    for st = 1:6
        
        betalim = [min(thesebetas(st,:)) max(thesebetas(st,:))];
        betalim = [-max(abs(betalim)) max(abs(betalim))];
        betarange = linspace(betalim(1),betalim(2),100);
        
        subplot(1,6,st)
        for i = 1:length(chanidx)
            scatter(pos(chanidx(i),1),pos(chanidx(i),2),25,'markerfacecolor',cmap(findMin(betarange,thesebetas(st,i)),:),'markeredgecolor','k'); hold on
        end
        axis equal
        xlim([-.5 .5])
        ylim([-.5 .5])
        set(gca,'ticklength',[0 0])
        axis off
        drawnow
        
    end
    sgtitle(subjects{s})
    
    % predicted probability
    cmap = colours(6,'viridis');
    figure
    set(gcf,'position',[1 -149 1920 973])
    for st1 = 1:6
        
        % load classifier for this state
        thisL = best_lambdas(s,trainTimes==besttime);
        thisB = squeeze(classifier.betas(st1,:,thisL));
        thisI = classifier.intercepts(st1,thisL);
        
        % apply classifier to different time points from data for different state trials
        pred = [];
        for t = 1:nT
            
            [X,Y] = organise_FL_data(data,trainTimes(t),1);
            
            % remove null data
            X = X(1:size(X)/2,:);
            Y = Y(1:size(X)/2,:);
            
            if trainTimes(t)==besttime % same data
                load(fullfile('D:\2020_RiskyReplay\data\meg\classifiers',subjects{s},['cv_' subjects{s} '_t' num2str(besttime) '_n1.mat']));
                for st2 = 1:6
                    pred(t,st2) = mean(squeeze(avpred(:,st1,st2,thisL)));
                end
            else            
                for st2 = 1:6
                    idx = Y == st2;
                    pred(t,st2) = mean(1 ./ (1 + exp(-(X(idx,:)*thisB' + repmat(thisI, [size(X(idx,:),1) 1])))));
                end
            end
        end
        
        % plot
        subplot(1,6,st1)
        for st2 = 1:6
            plot(trainTimes,pred(:,st2),'color',cmap(st2,:),'linewidth',1.4); hold on
        end
        plot(trainTimes([1 end]),[1/6 1/6],'k:'); hold on
        ylim([0 1])
        set(gca,'ticklength',[0 0])
        
    end
    sgtitle(subjects{s})
    
end

%% Cross-validate across time points

onCluster = true;

for s = 1:N
    
    % Load data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_FL_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    data = merged;
    clear merged
    
    % Remove incorrect trials
    behav = readtable(fullfile(dir_behav,subjects{s},[num2str(str2double(subjects{s})) '_fl.csv']));
    idx = strcmp(behav.Acc,'True'); % index of trials to include
    
    % Remove outlier RT trials
    idx = idx & abs(zscore(behav.RT)) < 5;
    
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

    % Create output directory for saving data/labels for the cluster
    dir_save = fullfile(dir_batch,subjects{s});
    if ~exist(dir_save)
        mkdir(dir_save)
    end
    
    % Cross-validate on other states
    for n = 1:nN
        cc = 0;
        for t1 = 1:nT
            for t2 = 1:nT
           
                if t1~=t2

                    disp('===================================')
                    disp(['=== SUBJECT ' subjects{s} ', TP1 = ' num2str(trainTimes(t1)) ' ms, TP2 = ' num2str(trainTimes(t2)) ' ms ==='])
                    disp('===================================')
    
                    % Load the data at T1 and T2
                    [X1,Y1] = organise_FL_data(data,trainTimes(t1),nulldata(n));
                    [X2,Y2] = organise_FL_data(data,trainTimes(t2),nulldata(n));
        
                    % Save 
                    cc = cc+1;
                    filename = fullfile(dir_save,['data_' subjects{s} '_idx' num2str(cc) '_n' num2str(nulldata(n)) '.mat']);
                    timeidx = trainTimes([t1 t2]);
                    save(filename,'X1','Y1','X2','Y2','timeidx');
    
                    if ~onCluster
                        crossvalidate_classifier(filename); % save to disk
                    end

                end
            end
        end

        if onCluster

            % Create job for cluster
            generate_jobs_classifier(subjects{s},[num2str(1) '-' num2str(cc)],nulldata(n),'holly');

        end

    end
end

%% After computing all the cross-validation, put into matrix
CV = nan(N,nT,nT,nStates);
for s = 1:N

    for t1 = 1:nT

        % Load classifier trained on T1
        load(fullfile('D:\2020_RiskyReplay\data\meg\classifiers',subjects{s},...
            ['classifier_' subjects{s} '_t' num2str(trainTimes(t1)) '_n' num2str(nulldata(n)) '.mat']));

        % load original accuracy for this classifier (same training/test time, so need to do folds)
        load(fullfile('D:\2020_RiskyReplay\data\meg\classifiers',subjects{s},...
            ['cv_' subjects{s} '_t' num2str(trainTimes(t1)) '_n' num2str(nulldata(n)) '.mat']));

        cv = squeeze(mean(cv));

        % get best lambda
        [~,bestLambda] = max(mean(cv),[],2);

        % log
        CV(s,t1,t1,:) = cv(:,bestLambda);

        for t2 = 1:nT

            if t1~=t2

                % load accuracy of t1 applied to t2 data


            end
        end
    end
end


plotcv = interp2(squeeze(mean(mean(CV),4)),5);

figure
imagesc(plotcv)
colormap(colours(256,'inferno'))
set(gca,'ticklength',[0 0])
