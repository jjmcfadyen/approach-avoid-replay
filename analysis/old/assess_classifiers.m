% Plot classification accuracy

clear all
clc

nType = 1;         % null data

%% Directories

addpath('utils')
[sinfo, dinfo] = dir_cfg();

addpath(dinfo.tb_fieldtrip);
ft_defaults

dir_data = fullfile(dinfo.data_meg_classifiers,['1vsRest_n' num2str(nType)]); % <-- choose iteration here

%% Parameters

subjects = sinfo.subjects;
N = length(subjects);

nL = 100; % no. of lambda values

nS = 6; % no. of states

imageNames = {'baby', 'backpack', 'bicycle', 'bowtie', 'car', 'cat',...
             'cupcake', 'hourglass', 'house', 'lamp', 'toothbrush', 'zebra'};
nI = length(imageNames);

%% Load data

% Set parameters
trainTimes = 0:10:300; % in ms
nT = length(trainTimes);

% Load data
predictions = nan(N,nT,nS,nS,nL); % subjects, training time, classifier, actual state, lambda --> prediction, 0 to 1
outcomes    = nan(N,nT,nS,nS,nL); % subjects, training time, classifier, actual state, lambda --> % of time winner across folds
accuracy    = nan(N,nT,nS,nL);    % subjects, training time, classifier, lambda --> accuracy across folds
for s = 1:N
   
    subject = subjects{s};
    
    for tb = 1:length(trainTimes)
        
        try
            load(fullfile(dir_data,subject,[subject '_cc_train-' num2str(trainTimes(tb)) 'ms.mat'])); % loads 'output' variable
        
            % the predictions of each classifier for each state
            predictions(s,tb,:,:,:) = squeeze(mean(output.predictions)); % average over folds

            % how accurate each classifer was on each fold
            accuracy(s,tb,:,:) = squeeze(mean(output.accuracy)); % average over folds

            % confusion matrix
            for classifier = 1:nS
                for st = 1:nS
                    outcomes(s,tb,classifier,st,:) =  mean(squeeze(output.outcomes(:,classifier,:)) == st);
                end
            end
        catch
            warning(['Unable to load file for ' subject ', ' num2str(trainTimes(tb)) ' ms']); 
        end
    end
end

%% HEATMAPS: Training Time x Lambda

figure
cmap = colours(length(0:.01:1),'inferno');
for s = 1:N+1
    
    subplot(4,4,s)
    
    if s <= N % individual subject
        
        % average over states
        Y = squeeze(mean(accuracy(s,:,:,:),3));
        title(subjects{s})
        
    else % grand average
        
        Y = squeeze(mean(mean(accuracy),3));
        title('Grand Average')
        colorbar
        
    end
    
    imagesc(Y');
    colormap(cmap);
    set(gca,'ticklength',[0 0])
    set(gca,'xtick',[0:5:nT])
    set(gca,'xticklabels',strsplit(num2str(trainTimes(1:5:nT))))
    ylabel('Lambda #')
    xlabel('Training Time (ms)')
    
    caxis([(1/6) .6])
    
end

%% LINE PLOT: Best lambda per subject --> images & subjects over training times

% --- get image names for each subject
imageCount = array2table(nan(N,length(imageNames)),'variablenames',imageNames); % 12 possible images
stateCount = nan(N,nS); % 6 states
stateNames = cell(N,nS);

for s = 1:N
    [fl,event,~] = pp_cfg(subjects{s},'FL');
    [training,names,~] = getFL(subjects{s},fl,event);
    for st = 1:nS
        imageCount(s,strfindcell(imageNames,names{st})) = array2table(sum(training.L == st));
        stateCount(s,st) = sum(training.L == st);
    end
    stateNames(s,:) = names;
end

% --- get best lambda per subject/time/classifier
bL = nan(N,nT,nI);
Y_img = nan(N,nT,nI);
Y_subj = nan(N,nT,nS);

for s = 1:N
    for t = 1:nT
        for i = 1:nS

            y = squeeze(accuracy(s,t,i,:));
            max_L = find(y == max(y));
            max_L = max_L(randi(length(max_L)));
            y = y(max_L);

            idx = strfindcell(imageNames,stateNames{s,i});
            Y_img(s,t,idx) = y;

            Y_subj(s,t,i) = y;
            bL(s,t,i) = max_L;
        end
    end
end

% --- plot
figure
for g = 1:2 % images, subjects

    subplot(1,2,g) % images

    if g == 1
        cmap = colours(nI,'viridis');
        m = squeeze(nanmean(Y_img));
        sem = squeeze(nanstd(Y_img));
        n = sum(~isnan(table2array(imageCount)));
        for i = 1:nI
            sem(:,i) = sem(:,i) / sqrt(n(i));
        end
    elseif g == 2
        cmap = colours(N,'plasma');
        m = squeeze(mean(Y_subj,3));
        sem = [];
        for s = 1:N
            sem(s,:) = std(squeeze(Y_subj(s,:,:))') / sqrt(nS);
        end
        m = m';
        sem = sem';
    end
    lower = m-sem;
    upper = m+sem;

    % reorder by max accuracy
    [~,sortidx] = sort(max(m),'ascend');
    m = m(:,sortidx);
    upper = upper(:,sortidx);
    lower = lower(:,sortidx);

    % plot individual
    for i = 1:size(m,2)
        plot(trainTimes,m(:,i),'color',cmap(i,:),'linewidth',1.4); hold on
        for t = 1:nT
            plot(repmat(trainTimes(t),2,1),[lower(t,i) upper(t,i)],'color',cmap(i,:),'linewidth',1.2,'handlevisibility','off'); hold on
        end
    end

    % plot grand average
    plot(trainTimes,mean(m,2),'k','linewidth',2.2); hold on

    % decorate
    title(['Max avg = ' num2str(round(max(mean(m,2))*100,2)) '% at ' num2str(trainTimes(find(mean(m,2) == max(mean(m,2))))) ' ms'])
    set(gca,'ticklength',[0 0])
    plot(trainTimes([1 end]),[1/6 1/6],'k--'); hold on

    if g == 1
        this_legend = imageNames(sortidx);
        for i = 1:length(this_legend)
            this_legend{i} = [this_legend{i} ' (' num2str(sum(~isnan(table2array(imageCount(:,sortidx(i)))))) ')'];
        end
        legend(this_legend)
    elseif g == 2
        this_legend = subjects(sortidx);
        for i = 1:length(this_legend)
            this_legend{i} = this_legend{i}(1:3);
        end
        legend(this_legend)
    end

end

%% CONFUSION MATRIX

% First, organise the data
        
bL = nan(N,nT,nS); % best lambdas per classifier

Y = nan(N,nT,nS,nS);
for s = 1:N
    for t = 1:nT
        for classifier = 1:nS
            y = squeeze(accuracy(s,t,classifier,:));
            max_y = find(y == max(y));
            max_y = max_y(randi(length(max_y)));
            bL(s,t,classifier) = max_y;
            Y(s,t,classifier,:) = squeeze(outcomes(s,t,classifier,:,max_y));
        end
    end
end

% concatenate into massive grid
G = [];
for s = 1:N+1

    g = [];
    for t = 1:nT
        if s <= N
            g = [g, squeeze(Y(s,t,:,:))];
        else
            g = [g, squeeze(mean(Y(:,t,:,:)))]; % average across subjects
        end
        if t < nT
            g = [g, zeros(nS,1)]; % add black space between grids
        end
    end

    G = [G; g];
    if s < N+1
       G = [G; zeros(1,size(G,2))]; 
    end

end

% plot
figure
cmap = colours(1000,'inferno');
imagesc(G)
colormap(cmap)
caxis([1/6 1])
set(gca,'visible','off')
    