%% Get temporal modulation betas

% Use the replay optimised training times from the planning period
load(fullfile('D:\2020_RiskyReplay\data\meg\replay\withoutintercept','optimised_times.mat'))

trialcount = nan(N,1);
avmodulation = [];
X = cell(1,N);
Y = cell(1,N);
for s = 1:N
   
    disp(['Getting replay for ' subjects{s} '...'])
    
    % Get replay from planing period
    load(fullfile('D:\2020_RiskyReplay\data\meg\replay\withoutintercept',...
        subjects{s},['replay_norm_time-modulated_' subjects{s} '_t' num2str(optimised_times(s)) '_n1.mat']));
    modulation = replay(:,:,4:6,:,:);
    
    % Log to group (averaged across trials & transitions)
    avmodulation(s,:,:,:) = squeeze(mean(mean(modulation(:,:,:,:,:),4)));
    
    % match up the replay data to the behavioural data
    load(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']))
    behav = behav.task;
    
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    trialinfo = merged.trialinfo;
    clear merged
    
    idx = nan(size(modulation,1),1); % converts the [block trial] index into the row index in 'behav'
    for trl = 1:length(idx)
        if trialinfo(trl,1)==0
            idx(trl,1) = find(behav.Practice==1 & behav.Trial==trialinfo(trl,2));
        else
            idx(trl,1) = find(behav.Practice==0 & behav.Block==trialinfo(trl,1) & behav.Trial==trialinfo(trl,2));
        end
    end
    
    % crop the 'behav' variable to only include those with an associated replay trial (e.g., some MEG blocks are missing)
    behav = behav(sort(idx),:);
    
    % sort replay trials in order of behavioural trials
    [~,sortidx] = sort(idx);
    modulation = modulation(sortidx,:,:,:,:);
    
    % Put in Y variable
    Y{s} = behav(:,[1:3 5:6 8:10 23:end-2]); % just take the most relevant columns, to save space
    
    %% Organise replay per path/choice/outcome
    
    % sort into good vs bad paths
    nTrls   = size(modulation,1); % trials
    nPerms  = size(modulation,2); % null permutations
    nSeq    = size(modulation,3); % forward, backward, forward-backward
    nTrans  = size(modulation,4); % number of transitions
    nLag    = size(modulation,5); % no. of lags (10:10:600 ms)
    
    rewarding = nan(nTrls,nPerms,nSeq,nLag); 
    aversive = nan(nTrls,nPerms,nSeq,nLag); 
    
    thisidx = behav.nV_1 > behav.nV_2 & behav.Forced==0; % trials where path 1 is better than path 2
    rewarding(thisidx,:,:,:) = squeeze(mean(modulation(thisidx,:,:,1:2,:),4));
    
    thisidx = behav.nV_1 < behav.nV_2 & behav.Forced==0; % trials where path 2 is better than path 1
    rewarding(thisidx,:,:,:) = squeeze(mean(modulation(thisidx,:,:,3:4,:),4));
    
    thisidx = behav.nV_1 < behav.nV_2 & behav.Forced==0; % trials where path 1 is worse than path 2
    aversive(thisidx,:,:,:) = squeeze(mean(modulation(thisidx,:,:,1:2,:),4));
    
    thisidx = behav.nV_1 > behav.nV_2 & behav.Forced==0; % trials where path 2 is worse than path 1
    aversive(thisidx,:,:,:) = squeeze(mean(modulation(thisidx,:,:,3:4,:),4));
    
    % Put in X variable
    x = nan(nTrls,2,nPerms,nSeq,nLag);
    x(:,1,:,:,:) = rewarding;
    x(:,2,:,:,:) = aversive;

    X{s} = x;
    
end

%% Plot group average, optimised training time

excludeSubjects = {'263098','680913'}; %{'396430'}; % bad classification
thesesubjects = find(~ismember(subjects,excludeSubjects));
thisN = length(thesesubjects);  

figure;
cc = 0;
for g = 1:3
    cc = cc+1;
    subplot(1,3,cc)
    opts = [];
    opts.g = 3;
    opts.nullThresh = false;
    opts.nullMoving = true;
    d = squeeze(avmodulation(thesesubjects,:,g,:));
    ss_plot(d,opts);
    if g==10
        title('forward');
    elseif g==2
        title('backward');
    elseif g==3
        title('fwd-bwd');
    end
end

%% Plot path-specific replay per choice/outcome

pathDefinition = 'pos-neg'; % 'pos-neg' (only trials where one path is better than the other) 
                            % 'better-worse' (all trials, including those where both are positive or both are negative)
cmap = [0, 238, 144
        255, 0, 89
        125 125 125]/255;

% Group average                            
pdata = [];
for s = 1:thisN
    
    B = Y{thesesubjects(s)};
    R = X{thesesubjects(s)};
    
    % remove forced-choice trials
    ridx = B.Forced==0;
    B = B(ridx,:);
    R = R(ridx,:,:,:,:);
    
    for c = 1:2 % choice
        idx = B.Choice==c;
        if strcmp(pathDefinition,'pos-neg')
            idx = idx & sum([B.nV_1 B.nV_2] > 1,2) == 1;
        end
        pdata(s,c,:,:,:,:) = squeeze(mean(R(idx,:,:,:,:)));
    end
end

figure
for i = 1:2
    opts = [];
    opts.avCol = cmap(i,:);
    opts.nullThresh = false;
    opts.nullMoving = true;
    cc = 0;
    for c = 1:2
        for g = [3 1 2]
            cc = cc+1;
            subplot(2,3,cc)
            opts.g = 3;
            d = squeeze(pdata(:,c,i,:,g,:));
            ss_plot(d,opts);
        end
    end
end

figure
cc = 0;
for c = 1:2
    opts = [];
    opts.avCol = cmap(3,:);
    opts.nullThresh = false;
    opts.nullMoving = true;
    for g = [3 1 2]
        cc = cc+1;
        subplot(2,3,cc)
        opts.g = 3;
        d = squeeze(pdata(:,c,1,:,g,:)) - squeeze(pdata(:,c,2,:,g,:));
        ss_plot(d,opts);
    end
end

%% Linear mixed effects modelling

subtractNull = false;
g = 3;

% Get replay
replay_window = 2:9; % lags to include
lags = 10:10:600; % in ms

T = [];
lmetrialcount = nan(thisN,1);
for s = 1:thisN
	
    % pick trials
    idx =  Y{s}.Forced==0;
    if strcmp(pathDefinition,'pos-neg')
        idx = idx & sum([Y{s}.nV_1 Y{s}.nV_2] > 1,2) == 1;
    end
    
    B = Y{s}(idx,:);
    
    % get replay
    if subtractNull
        d = [];
        for i = 1:2
            tmp = squeeze(X{s}(idx,i,:,g,:));
            np = squeeze(mean(tmp(:,2:end,:)));
            if g==3
                abs_np = unique(round(abs(np),4),'rows');
                npthresh = quantile(max(abs_np,[],2),.975);
                d(:,i,:) = abs(squeeze(tmp(:,1,:))) > npthresh;
                posidx = squeeze(tmp(:,1,:)) > 0;
                for j = 1:size(d,1)
                    d(j,i,~posidx(j,:)) = d(j,i,~posidx(j,:)) * -1;
                end
            else
                np = unique(round(np,4),'rows');
                npthresh = quantile(max(np,[],2),.95);
                d(:,i,:) = squeeze(tmp(:,1,:)) > npthresh;
            end
        end
    else
        d = squeeze(X{s}(idx,:,1,g,replay_window)); % get just this subject at the specific lags
    end

    d_rewarding = squeeze(d(:,1,:));
    d_aversive = squeeze(d(:,2,:));
    d_differential = d_rewarding - d_aversive;
    
    lmetrialcount(s,1) = size(B,1);
    
    % combine with behavioural data
    n = size(d,1);
    thistable = [];
    for t = 1:length(replay_window)
        tmp = B;
        tmp.Lag = repmat(lags(replay_window(t)),n,1);
        tmp.Modulation_rewarding = d_rewarding(:,t);
        tmp.Modulation_aversive = d_aversive(:,t);
        tmp.Modulation_differential = d_differential(:,t);
        if t==1
            thistable = tmp;
        else
            thistable = [thistable; tmp];
        end
    end
    
    thistable.Subject = repmat(subjects(thesesubjects(s)),size(thistable,1),1);
    
    if s>1
        thistable = thistable(:,contains(thistable.Properties.VariableNames,T.Properties.VariableNames));
    end
    
    T = [T; thistable];
    
end

% Recode variables
T.Choice = T.Choice+10;
T.Choice(T.Choice==11) = 1; % approach
T.Choice(T.Choice==12) = 0; % avoid

T.Acc = categorical(T.Acc,[0 1],{'incorrect','correct'});
T.Subject = categorical(T.Subject,unique(T.Subject),unique(T.Subject));

% Factors that predict choice
lme = fitlme(T(T.Choice==1,:),'Modulation_aversive~EV+(1|Subject:Lag)')

% Save
writetable(T,'D:\2020_RiskyReplay\results\replay\replay_temporal_lme.csv');
