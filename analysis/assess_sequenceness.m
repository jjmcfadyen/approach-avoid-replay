% Assess sequenceness

clear all
clc

e = 1; % 1 = decision, 2 = transition, 3 = image, 4 = outcome

decisionType = 'full'; % 'tails', or 'full' decision epoch
decisionRange = 'half'; % in seconds OR 'half'

startTime = 0; % in seconds, where to take the epoch from

%% Directories

addpath('utils')
[sinfo, dinfo] = dir_cfg();

addpath(dinfo.tb_fieldtrip);
ft_defaults;

dir_seq = fullfile(dinfo.results,'meg','sequenceness');
dir_load = fullfile(dinfo.data_meg,'seqOptimisation','1vsRest_n1');

%% Parameters

subjects = sinfo.subjects;
N = length(subjects);

trainTimes  = 100:10:300; % in ms
timerange   = 120:10:150; % time windows used for optimisation
useTime     = [120 130 140]; % which training time(s) to use (if multiple, they are averaged)

epochs = {'decision','transition','image','outcome'};
nE = length(epochs);

lags = 10:10:600;

% Get all relevant null permutations
opts = [];
opts.pathLink = false;
opts.withinSeq = false;
U = ss_nullPerms(opts); 

% Only keep unique transitions, ignoring order
opts = [];
opts.ignoreOrder = true;
opts.removePos = false;
opts.removeLoop = false;
uIdx = ss_editNull(U,opts);
rU = U(uIdx,:); % this is what was done with the sequence optimisation

epochs = {'decision','transition','image','outcome'};
epoch = epochs{e};

%% Rearrange data
   
load(fullfile(dir_seq,['seq_bestLambdas_decision.mat'])); % loads 'cBL' and 'CA' variables

% get sequenceness per subject, creating trial-stacked array
sq = [];
behav = [];
spred = cell(1,N); % save predicted data
for s = 1:N
    
    subject = subjects{s};
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf(['Loading sequenceness data for ' subject '...\n'])
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    
    t = find(trainTimes == timerange(stimewin(s)));
    this_lambda = cBL(s,t);
    
    % load info
    [fl,event,in] = pp_cfg(subject,'Task');
    [test, ~] = getTask(subject,fl,event,epoch);
        
    % Build classifier & extract optimal lambda
    classifier = cc_build(subject,trainTimes(t),1,'1vsRest');
    classifier.B = classifier.B(:,:,this_lambda);
    classifier.I = classifier.I(:,this_lambda)';

    % Load task data for this subject & epoch
    in.Fs = test.Fs;
    nTrls = size(test.L,1);

    % Get behaviour
    thisT = event.table;
    if e == 1
        idx = ismember([thisT.ExpTrial thisT.Choice],test.L,'rows');
    elseif e == 2
        idx = ismember([thisT.ExpTrial thisT.Transition],test.L,'rows');
    elseif e == 4
        idx = ismember(thisT.ExpTrial,test.L(:,1),'rows');
    end
    if sum(idx) ~= nTrls
        error([subject ' trials did not match between behaviour and MEG'])
        missing = ~ismember(test.L,[thisT.ExpTrial thisE],'rows');
        test.L(missing,:) = [];
        test.x(:,missing) = [];
        test.D(:,missing) = [];
        nTrls = size(test.L,1);
    end
    thisT = thisT(idx,:);
    thisT.Subject = repmat({subject},size(thisT,1),1);
    if s == 1
        behav = thisT; 
    else
        behav = [behav; thisT];
    end

    % Edit start time (if applicable)
    if startTime > 0
        for trl = 1:nTrls
             test.D{trl} = test.D{trl}(test.x{trl} > startTime,:);
             test.x{trl} = test.x{trl}(test.x{trl} > startTime);
        end
    end

    % Generate sequenceness
    this_sq = [];
    spred{s} = cell(1,nTrls);
    for trl = 1:nTrls

        if strcmp(epoch,'decision') && strcmp(decisionType,'tails') % take tails 
            for i = 1:2 % first part, second part
                
                if ischar(decisionRange)
                    midpoint = round(linspace(1,size(test.D{trl},1),3));
                    if i ==1
                        elim = [1 midpoint(2)-1]; % start
                    elseif i == 2
                        elim = [midpoint(2) midpoint(end)]; % end
                    end
                else
                    if i ==1
                        elim = [1 find(test.x{trl} < decisionRange,1,'last')]; % start
                    elseif i == 2
                        elim = [find(test.x{trl} >= test.x{trl}(end)-decisionRange,1,'first') length(test.x{trl})]; % end
                    end
                end

                thisD = test.D{trl}(elim(1):elim(2),:);
                pred = cc_predict(thisD,classifier,1,1);
                this_sq(trl,:,:,:,:,i) = ss_build(pred,U,in);
                spred{s}{i,trl} = pred;
            end
        else % generate sequenceness for other epoch (transition, image, outcome)
            pred = cc_predict(test.D{trl},classifier,1,1);
            this_sq(trl,:,:,:,:) = ss_build(pred,U,in);
            spred{s}{trl} = pred;
        end            
    end
    
    if s == 1
        sq = this_sq;
    else
        t0 = size(sq,1)+1;
        t1 = size(this_sq,1)+(t0-1);
        if e==1 && strcmp(decisionType,'tails')
            sq(t0:t1,:,:,:,:,:) = this_sq;
        else
            sq(t0:t1,:,:,:,:) = this_sq;
        end
    end
end
clear this_sq
clear tmp
clear pred

if strcmp(epoch,'decision') && strcmp(decisionType,'tails')
    suffix = '_half';
else
    suffix = '';
end

save(fullfile(dir_seq,['behavTable_' epochs{e} '.mat']),'behav');
save(fullfile(dir_seq,['sq_' epochs{e} suffix '.mat']),'sq');
save(fullfile(dir_seq,['pred_' epochs{e} suffix '.mat']),'spred');

%{
load(fullfile(dir_seq,['behavTable_' epochs{e} '.mat']));
load(fullfile(dir_seq,['sq_' epochs{e} suffix '.mat']));
load(fullfile(dir_seq,['pred_' epochs{e} suffix '.mat']));
%}

%% Overall sequenceness

T = behav;
D = sq;

% fix permutations (the sequence optimisation was done on only unique
% transitions, regardless of order, but now we want to compare specific
% transitions with different orders)
tp = [1 2; 2 3; 4 5; 5 6];
cp = [];
for i = 1:size(U,1)
    thisPerm = U(i,:);
    computed_perm = find(ismember(rU,thisPerm,'rows'));
    if ~isempty(computed_perm)
        cp = [cp; computed_perm];
    else
        % within the reduced permutations (used for optimisation), find the
        % permutation that matches this full permutation (just in a
        % different order)
        computed_perm = find(ismember(rU(:,1:3),U(i,4:6),'rows') & ismember(rU(:,4:6),U(i,1:3),'rows'));
        matching_perm = rU(computed_perm,:);
        for t = 1:size(tp,1)
            idx = [find(thisPerm==matching_perm(tp(t,1))) find(thisPerm==matching_perm(tp(t,2)))];
            idx = all(ismember(tp,idx),2);
            D(:,i,:,t,:,:) = squeeze(D(:,computed_perm,:,idx,:,:));
        end
        cp = [cp; computed_perm];
    end
end

% average across subjects
sidx = setdiff(1:N,[]);
n = length(sidx);
lagRange = lags <= 600;

y = [];
for s = 1:n
    idx = strfindcell(T.Subject,subjects{sidx(s)}) & T.Forced==0;
    y(s,:,:,:,:,:) = squeeze(mean(mean(D(idx,:,:,:,lagRange,:),1),4));
end

opts = [];
opts.showIndividual = false;
opts.indCol = colours(size(y,1),'plasma');
opts.x = lags(lagRange);
figure
cc = [1:3;4:6];
for g = 1:3
    opts.g = g;
    if e==1 && strcmp(decisionType,'tails')
        for i = 1:2
            subplot(2,3,cc(i,g))
            ss_plot(squeeze(y(:,:,g,:,i)),opts);
        end
    else
        subplot(1,3,g)
        ss_plot(squeeze(y(:,:,g,:)),opts);
    end
end

%% Get sequenceness for each path

choosePType = 'relative'; % 'absolute' 
                          %     (only select trials where the rewarding path 
                          %     was positive and the averisve path was negative, but lose 
                          %     all catch trials ~ 15% of trials), or...
                          % 'relative' <--- THIS IS WHAT WAS USED IN THE PAPER
                          %     (denote a path as 'rewarding' if it was
                          %     bigger than the other path, or 'aversive'
                          %     if it was smaller than the other path, even
                          %     if both are positive or both are negative)

nTrls = size(D,1);

% identfy path types
ptype = nan(nTrls,2); % 1st col = path 1 type, 2nd col = path 2 type (1 = rewarding, -1 = aversive)

if strcmp(choosePType,'absolute')
    
    ptype(T.nV_1 > 0,1) = 1;  % path 1 rewarding
    ptype(T.nV_1 <= 0,1) = -1; % path 1 aversive
    ptype(T.nV_2 > 0,2) = 1;  % path 2 rewarding
    ptype(T.nV_2 <= 0,2) = -1; % path 2 aversive

    % remove trials where BOTH were aversive or BOTH were rewarding
    idx = ptype(:,1)~=ptype(:,2);
    ptype = ptype(idx,:);
    T = T(idx,:);
    D = D(idx,:,:,:,:,:);

elseif strcmp(choosePType,'relative')
    
    ptype(T.nV_1 > T.nV_2,1) = 1; % path 1 more rewarding, path 2 more aversive
    ptype(T.nV_1 > T.nV_2,2) = -1; 
    ptype(T.nV_1 < T.nV_2,1) = -1; % path 2 more rewarding, path 1 more aversive
    ptype(T.nV_1 < T.nV_2,2) = 1; 
    
    idx = T.nV_1 ~= T.nV_2; % sometimes the values of each path are exactly the same
    T = T(idx,:);
    D = D(idx,:,:,:,:,:);
    ptype = ptype(idx,:);
    
end

nTrls = size(D,1);

% create new sequenceness variable for each path, updating the null as well
ds = nan(nTrls,size(D,2),size(D,3),2,size(D,5),size(D,6)); % trials, perms, g, path, lag
actualPath = nan(nTrls,1); % col1 = rewarding, col2 = aversive
for trl = 1:nTrls
    
    % extract sequenceness for each permutation
    y = nan(size(D,2),3,2,size(D,5),size(D,6)); % perms, g, path type, lags
    for a = 1:2
        
        % get transitions for this path (a = actual physical path number)
        if a == 1
            idx = [1 2];
        elseif a == 2
            idx = [3 4];
        end
        
        if ptype(trl,a) == 1 % if this path is rewarding
            y(:,:,1,:,:) = squeeze(mean(D(trl,:,:,idx,:,:),4)); % average across path transitions
            actualPath(trl,1) = a;
        elseif ptype(trl,a) == -1 % if this path is aversive
            y(:,:,2,:,:) = squeeze(mean(D(trl,:,:,idx,:,:),4)); % average across path transitions
            actualPath(trl,2) = a;
        end
        
    end
    
    % insert into variable
    ds(trl,:,:,:,:,:) = y;
 
end

% Plot to check (should be the same as before)
sidx = setdiff(1:N,[]);
n = length(sidx);
lagRange = lags <= 600;

y = [];
for s = 1:n
    idx = strfindcell(T.Subject,subjects{sidx(s)}) & T.Forced==0;
    y(s,:,:,:,:) = squeeze(mean(mean(ds(idx,:,:,:,lagRange,:),1),4));
end

figure
opts = [];
opts.x = lags(lagRange);
opts.showIndividual = false;
opts.indCol = colours(size(y,1),'plasma');
cc = [1:3; 4:6];
for g = 1:3
    opts.g = g;
    if e==1 && strcmp(decisionType,'tails')
        for i = 1:2
            subplot(2,3,cc(i,g))
            ss_plot(squeeze(y(:,:,g,:,i)),opts);
        end
    else
        subplot(1,3,g)
        ss_plot(squeeze(y(:,:,g,:)),opts);
    end
end

%% Plot sequenceness for rewarding/aversive paths

lagRange = lags <= 300;

% select subjects
sidx = setdiff(1:N,[2 3]);
n = length(sidx);

if e==1 && strcmp(decisionType,'tails')
    I = 2;
    ttext = {'Start decision: ','End decision: '};
else
    I = 1;
    ttext = {epochs{e}};
end

% plot rewarding vs aversive paths per choice
ysubj = cell(n,2);
for i = 1:I
    
    tC = [];
    for c = 1:2 % decision: choices, transition: outcome
        
        figure
        if e == 1 || e == 4
            if c == 1
                sgtitle([ttext{i} 'Risky'])
            elseif c == 2
                sgtitle([ttext{i} ': Safe'])
            end
        elseif e == 2
            if c == 1
                sgtitle([ttext{i} ': Selected'])
            elseif c == 2
                sgtitle([ttext{i} ': Unselected'])
            end
        end

        y = [];
        for s = 1:n
            subject = subjects{sidx(s)};
            if e == 1 || e == 4 % get rewarding & aversive paths for this choice
                idx = strfindcell(T.Subject,subject,'logical') & T.Choice == c & T.Forced==0;
                y(s,:,:,:,:) = squeeze(mean(ds(idx,:,:,:,lagRange,i)));
                ysubj{s,c} = squeeze(ds(idx,:,:,:,lagRange,i));
            elseif e == 2 % get rewarding and aversive paths for this type of transition (selected vs unselected)
                idx = strfindcell(T.Subject,subject,'find');
                tmp = cell(1,2); % rewarding, aversive
                cc = [0 0];
                for trl = 1:length(idx)
                    thisT = T(idx(trl),:); 
                    selected = find(actualPath(idx(trl),:) == thisT.Transition);
                    unselected = setdiff(1:2,selected);
                    if c == 1 % selected
                        cc(selected) = cc(selected) + 1;
                        tmp{selected}(cc(selected),:,:,:) =  squeeze(ds(idx(trl),:,:,selected,:));
                    elseif c == 2 % unselected
                        cc(unselected) = cc(unselected) + 1;
                        tmp{unselected}(cc(unselected),:,:,:) =  squeeze(ds(idx(trl),:,:,unselected,:));
                    end
                end
                for k = 1:2
                    tC(s,c,k) = size(tmp{k},1);
                    y(s,:,:,k,:) = squeeze(mean(tmp{k}(:,:,:,lagRange)));
                end
            end
        end

        opts = [];
        opts.x = lags(lagRange);
        opts.twosided = true;
        opts.nullThresh = false;
        cmap = [18, 220, 161; 255, 49, 99; 138, 97, 255]/255; % rewarding, aversive, difference
        for p = 1:2
            opts.avCol = cmap(p,:);
            for g = 1:3
                subplot(2,3,g)
                opts.g = g;
                ss_plot(squeeze(y(:,:,g,p,:)),opts);
            end
        end

        opts.avCol = cmap(3,:);
        opts.twosided = true;
        for g = 1:3
            subplot(2,3,g+3)
            opts.g = g;
            ss_plot(squeeze(y(:,:,g,1,:)-y(:,:,g,2,:)),opts);
        end
    end
end


%% Plot by VALUE

V = [T.nV_1 T.nV_2];

if strcmp(choosePType,'absolute')
    
    rV = V(:,1);
    rV(rV < 0) = NaN;
    rV(isnan(rV)) = V(isnan(rV),2);

    lV = V(:,1);
    lV(lV > 0) = NaN;
    lV(isnan(lV)) = V(isnan(lV),2);
    
elseif strcmp(choosePType,'relative')
   
    rV = nan(size(V,1),1);
    lV = rV;
    
    [~,I] = max(ptype,[],2);
    I(:,2) = I(:,1) - 1;
    I(I(:,2) == 0,2) = 2;
    
    for i = 1:size(V,1)
        
        rV(i,1) = V(i,I(i,1));
        lV(i,1) = V(i,I(i,2));
        
    end
    
end

%% Plot sequenceness for different reward probabilities paths, across path types

differenceWaves = false;
collapseProbs = false; % turn 5 categories into 3 (unlikely, 50/50, likely)

% select subjects
sidx = setdiff(1:N,[2 3]);
n = length(sidx);

% calculate probability of rewarding & aversive paths
probs = nan(nTrls,2); % col1 = rewarding, col2 = aversive

idx = ptype(:,1) == 1; % when path 1 is rewarding...
probs(idx,1) = T.P(idx); % ... rewarding probability is P (path 1)
probs(idx,2) = 1-T.P(idx); % ... aversive probability is 1-P (path 2)

idx = ptype(:,1) == -1; % when path 1 is aversive...
probs(idx,1) = 1-T.P(idx); % ... rewarding probability is 1-P (path 2)
probs(idx,2) = T.P(idx); % ... aversive probability is P (path 1)

probs = round(probs,1);

% convert to categories
if collapseProbs
    probs(probs(:)==.1) = .3;
    probs(probs(:)==.9) = .7;
end

plist = unique(probs,'rows');

% check if comparing start/end of decision trials
if e==1 && strcmp(decisionType,'tails')
    I = 2;
    ttext = {'Start decision','End decision'};
else
    I = 1;
    ttext = {epochs{e}};
end

% plot
for i = 1:I
    
    tC = [];

    figure
    sgtitle(ttext{i})

    y = nan(n,length(plist),size(ds,2),size(ds,3),size(ds,4),size(ds,5));
    for s = 1:n
        for pr = 1:length(plist)
            idx = strfindcell(T.Subject,subjects{sidx(s)},'logical') & ismember(probs,plist(pr,:),'rows');
            tC(s,pr) = sum(idx);
            if any(idx)
                if sum(idx) == 1
                    y(s,pr,:,:,:,:) = ds(idx,:,:,:,:,i);
                else
                    y(s,pr,:,:,:,:) = squeeze(mean(ds(idx,:,:,:,:,i)));
                end
            else
                y(s,pr,:,:,:,:) = nan(size(ds,2),size(ds,3),size(ds,4),size(ds,5));
            end
        end
    end

    opts = [];
    opts.x = lags;
    opts.nullThresh = false;
    cmap = colours(length(plist),'viridis');
    for path = 1:2
        if ~differenceWaves || path == 1
            for pr = 1:length(plist)
                opts.avCol = cmap(pr,:);
                for g = 1:3
                    if ~differenceWaves
                        if path == 1
                            subplot(2,3,g)
                        elseif path == 2
                            subplot(2,3,g+3)
                        end
                        thisy = squeeze(y(:,pr,:,g,path,:));
                    else
                        subplot(1,3,g)
                        thisy = squeeze(y(:,pr,:,g,1,:)-y(:,pr,:,g,2,:));
                    end
                    opts.g = g;
                    ss_plot(thisy,opts);
                    if ~differenceWaves
                        if path == 1
                            ylabel('Rewarding sequenceness')
                        elseif path == 2
                            ylabel('Aversive sequenceness')
                        end
                    else
                        ylabel('aversive <--- ---> rewarding')
                    end
                end
            end
        end
    end
    L = legend(num2str(plist(:,1)));
    L.Title.String = 'Reward probability';
    
    % plot per probability condition
    cmap = [18, 220, 161; 255, 49, 99; 138, 97, 255]/255; % rewarding, aversive, difference
    opts = [];
    opts.x=lags;
    cc = [1 2 3; 4 5 6; 7 8 9];
    figure
    for pr = 1:length(plist)
        for g = 1:3
            opts.g = g;
            if ~differenceWaves
                for path = 1:2
                    subplot(3,3,cc(pr,g))
                    thisy = squeeze(y(:,pr,:,g,path,:));
                    opts.avCol = cmap(path,:);
                    ss_plot(thisy,opts);
                end
            else
                subplot(1,3,g)
                thisy = squeeze(y(:,pr,:,g,1,:)-y(:,pr,:,g,2,:));
                ss_plot(thisy,opts);
            end
            ylabel([num2str(plist(pr,1)*100) '%'])
        end
    end
        
end

%% Put into long format

nTrls = size(T,1);

% check if comparing start/end of decision trials
if e==1 && strcmp(decisionType,'tails')
    I = 2;
    suffix = {['-first' num2str(decisionRange)],['-last' num2str(decisionRange)]};
else
    I = 1;
    suffix = {''};
end

vN = [T.Properties.VariableNames 'PathProb' 'Direction' 'Path' 'Val' 'RewVal' 'LossVal' 'Seq' 'AboveThresh' 'Lag' 'AvSeq' 'AvSeqAvLag'];

if e == 2
    vN = [vN, 'Selected', 'Expected']; 
end

for i = 1:I
    avd = array2table(nan(0,length(vN)),'variablenames',vN);
    for g = 1:3
        for p = 1:3 % path type (3 = differential: rewarding - aversive, per trial)

            if p ==3
                seqdata = squeeze(ds(:,:,g,1,:,i)-ds(:,:,g,2,:,i));
            else
                seqdata = squeeze(ds(:,:,g,p,:,i));
            end

            avseqdata = squeeze(mean(ds,4));
            
            [~, npThresh] = so_subtractNull(seqdata,g,'trial');

            for l = 1:length(lags)

                tmp = T;

                tmp.Direction = repmat(g,nTrls,1);
                tmp.Lag = repmat(lags(l),nTrls,1);

                tmp.Path = repmat(p,nTrls,1);
                tmp.Seq = seqdata(:,1,l);
                tmp.AboveThresh = squeeze(seqdata(:,1,l)) - npThresh;
                
                tmp.AvSeq = squeeze(avseqdata(:,1,g,l));
                tmp.AvSeqAvLag = squeeze(mean(avseqdata(:,1,g,:),4));

                if p == 3
                    tmp.PathProb = probs(:,1); % for differential sequenceness, just give REWARDING probability
                    tmp.Val = rV-lV;
                else
                    tmp.PathProb = probs(:,p);
                    if p == 1
                        tmp.Val = rV;
                    elseif p == 2
                        tmp.Val = lV;
                    end
                end

                tmp.RewVal = rV;
                tmp.LossVal = lV;

                if e == 2
                    tmp.Selected = nan(size(tmp,1),1);
                    tmp.Expected = nan(size(tmp,1),1);
                    for trl = 1:size(tmp,1)
                        if p == 3
                            tmp.Selected(trl) = find(actualPath(trl,:) == tmp.Transition(trl)); % 1 = rewarding, 2 = aversive
                        else
                            tmp.Selected(trl) = double(actualPath(trl,p) == tmp.Transition(trl)); % true (1) or false (0)
                        end
                        if T.Transition(trl) == 1
                            tmp.Expected(trl) = T.P(trl);
                        elseif T.Transition(trl) == 2
                            tmp.Expected(trl) = 1-T.P(trl);
                        end
                    end
                end

                avd = [avd; tmp];
            end
        end
    end

    filename = fullfile(dir_seq,['longSeq_' epochs{e} suffix{i} '_withCatch.csv']);
    writetable(avd,filename);
    
end
