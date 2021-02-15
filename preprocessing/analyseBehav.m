%% analyse_behav.m
% Uses: data downloaded from Firestore via exportData Python script (on Mac)
%       - [SUBJECT]_s2_data.json
%       - [SUBJECT]_s2_structure.js
% Produces: data in table form per subject, plus group summary table and trial-by-trial long table:
%       - [SUBJECT]_task.csv
%       - groupBehaviour.mat
%       - longData.csv

clear all
close all
clc

addpath('utils');

%% Settings

plotIndividual = false;

%% Set variables

[sinfo, dinfo] = dir_cfg();

data = [];
data.structure = sinfo.sessions;
data.questionnaires = sinfo.questionnaires;

subjects = sinfo.subjects;
nC = sinfo.nC;

dir_data = dinfo.data_behav;
dir_results = dinfo.results_behav;

%% Run

T = [];
for subj = 1:length(subjects)
    
    subject = subjects{subj};

    fbehav = fullfile(dir_data,subject,[subject '_s2_data.json']);
    fstruct = fullfile(dir_data,subject,[subject '_s2_structure.js']);
    fsave = fullfile(dir_data,subject,[subject '_task.csv']);

    %% Import data

    % Behavioural output
    fid = fopen(fbehav); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    behav = jsondecode(str);

    % Experiment structure
    fid = fopen(fstruct); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    if strcmp(str(1),'v')
        str = str(16:end-2);
    end
    bstruct = jsondecode(str);

    %% Parameters

    ntrials = length(struct2array(bstruct.Trial));

    colnames = fieldnames(bstruct);
    binfo = array2table(nan(ntrials,length(colnames)),'VariableNames',colnames);
    for col = 1:5
        thisCol = extractfield(bstruct,colnames{col});
        x = struct2array(thisCol{1})';
        if size(x,2) > 1
            x = num2cell(x,2);
        else
            x = array2table(x);
        end
        binfo(:,col) = x; 
    end

    %% Analyse

    rNames = {'Practice','Block','Trial','Catch','Forced','ExpTrial','nCombo','P','nV_1','nV_2',...
              'S1','S2','S3','S4','S5','S6','nS1','nS2','nS3','nS4','nS5','nS6',...
              'EV','Choice','RT','Acc','Transition','Outcome'};
    results = array2table(nan(ntrials,length(rNames)),'VariableNames',rNames);
    for trl = 1:ntrials
        
        % get basic information (trial num, nCombo, P, nV_1, nV_2, EV)
        results(trl,1:5) = binfo(trl,1:5);
        results.ExpTrial(trl) = trl;

        nlist = struct2array(bstruct.N)';
        results.nCombo(trl) = find(nC(:,1) == nlist(trl,1) & nC(:,2) == nlist(trl,2));

        plist = struct2array(bstruct.P)';
        results.P(trl) = plist(trl,1);

        nVlist = reshape(extractfield(bstruct.nV,['x' num2str(trl-1)]),2,3);
        results.nV_1(trl) = nVlist(1,end);
        results.nV_2(trl) = nVlist(2,end);

        evlist = struct2array(bstruct.EV)';
        results.EV(trl) = evlist(trl,1);

        % see if is practice or not
        isPractice = results.Practice(trl)==1;
        
        % get path values (normal and flipped)
        Vlist = reshape(extractfield(bstruct.V,['x' num2str(trl-1)]),2,3);
        results(trl,11:16) = array2table([Vlist(1,:) Vlist(2,:)]);
        results(trl,17:22) = array2table([nVlist(1,:) nVlist(2,:)]);
        
        % get subject responses
        if isPractice
            tmp = behav.x5_negator_learning;
        else
            tmp = behav.x6_test;
        end

        idx = find(tmp.trial == extractfield(bstruct.Trial,['x' num2str(trl-1)]) & tmp.block == extractfield(bstruct.Block,['x' num2str(trl-1)]));
        if isempty(idx)
            warning(['Subject ' num2str(subject) ': could not find trial ' num2str(extractfield(bstruct.Trial,['x' num2str(trl-1)])) ' in block ' ...
                num2str(extractfield(bstruct.Block,['x' num2str(trl-1)]))])
            results(trl,:) = array2table(nan(1,size(results,2)));
        else
        
            if ~isempty(strmatch(tmp.choice{idx},'Airlock'))
                results.Choice(trl) = 1;
            elseif ~isempty(strmatch(tmp.choice{idx},'SUPPLY ROOM'))
                results.Choice(trl) = 2;
            end

            results.RT(trl) = tmp.rt(idx)/1000; % in seconds
            results.Acc(trl) = (results.Choice(trl) == 1 && results.EV(trl) >= 1) || (results.Choice(trl) == 2 && results.EV(trl) <= 1);
            if extractfield(bstruct.Forced,['x' num2str(trl-1)]) ~= 0
                results.Acc(trl) = NaN;
            end

            if results.Choice(trl) == 2
                results.Outcome(trl) = 1;
            else
                if ~isempty(strmatch(tmp.transition{idx},'DOOR 1'))
                    results.Outcome(trl) = results.nV_1(trl);
                elseif ~isempty(strmatch(tmp.transition{idx},'DOOR 2'))
                    results.Outcome(trl) = results.nV_2(trl);
                end
            end

            if strcmp(tmp.transition(idx),'DOOR 1')
                results.Transition(trl) = 1;
            elseif strcmp(tmp.transition(idx),'DOOR 2')
                results.Transition(trl) = 2;
            elseif strcmp(tmp.transition(idx),'SUPPLY ROOM')
                results.Transition(trl) = 0;
            end
        end
    end

    % remove all NaN rows
    idx = [];
    for trl = 1:size(results,1)
        tmp = isnan(table2array(results(trl,:)));
        if sum(tmp) == length(tmp)
            idx = [idx; trl];
        end
    end
    results(idx,:) = [];
    
    % Clean up other factors
    results.bAcc = results.Acc; % put 1 point buffer either side of EV to score accuracy
    results.bAcc(results.EV >= 0 & results.Choice == 1) = 1;
    results.bAcc(results.EV <= 2 & results.Choice == 2) = 1;
    results.bAcc(isnan(results.Acc)) = NaN;
    results.RT = results.RT + 5; % the first 5 seconds aren't included in the RT
    
    % Count up the trials that will be included
    includeIdx = results.Forced == 0 & results.RT < 30;
    results.Include = includeIdx;
    
    disp(['Subject ' num2str(subject) ': ' num2str(size(results,1)) ' trials in total'])
    
    writetable(results,fsave);
    data.tables{subj} = results;
    data.acc(subj,1) = nanmean(results.Acc(includeIdx));
    data.bAcc(subj,1) = nanmean(results.bAcc(includeIdx));
    
    thisRT = {};
    thisRT{1} = results.RT(includeIdx & results.Choice == 1);
    thisRT{2} = results.RT(includeIdx & results.Choice == 2);
    
    data.rt_mean(subj,:) = [mean(thisRT{1}) mean(thisRT{2})];
    data.rt_median(subj,:) = [median(thisRT{1}) median(thisRT{2})];
    data.rt_min(subj,:) = [min(thisRT{1}) min(thisRT{2})];
    data.rt_max(subj,:) = [max(thisRT{1}) max(thisRT{2})];
    
    % get indifference point
    [x, sortidx] = sort(results.EV(includeIdx));
    y = results.Choice(includeIdx)-1; % binary (0, 1) for choices (risky, safe)
    y = y(sortidx);
    [B,dev,stats] = glmfit(x,y,'binomial','logit');
    mfit = glmval(B,x,'logit');
    IP = mfit(find(abs(x - 1) == min(abs(x - 1)))); % indifference point
    data.ip(subj,1) = IP(1);

    tmp = results;
    tmp.Subject = repmat({subject},size(results,1),1);
    
    tmp.IUS = repmat(data.questionnaires(subj,1),size(tmp,1),1);
    tmp.Worry = repmat(data.questionnaires(subj,2),size(tmp,1),1);
    tmp.Ethical = repmat(data.questionnaires(subj,3),size(tmp,1),1);
    tmp.Financial = repmat(data.questionnaires(subj,4),size(tmp,1),1);
    tmp.Health = repmat(data.questionnaires(subj,5),size(tmp,1),1);
    tmp.Recreational = repmat(data.questionnaires(subj,6),size(tmp,1),1);
    tmp.Social = repmat(data.questionnaires(subj,7),size(tmp,1),1);
    
    T = [T; tmp];
    
    %% Plot

    if plotIndividual

        close all
        
        % %%%%%% 
        % %%%%%% Path outcomes and choices
        H1 = figure;

        % Path outcomes (dotted)
        % plot(results.ExpTrial,results.nV_1,'Color',[51 218 255]/255,'LineWidth',1.25,'LineStyle',':'); % path 1
        % hold on
        % plot(results.ExpTrial,results.nV_2,'Color',[255 141 51]/255,'LineWidth',1.25,'LineStyle',':'); % path 2
        ylim([min(results.Outcome)-1 max(results.Outcome)+1])

        % EV (solid)
        plot(results.ExpTrial,results.EV,'Color',[0 146 240]/255,'LineWidth',1.5,'LineStyle',':'); hold on

        % Participant choices
        ax = gca;

        x = results.ExpTrial(includeIdx & results.Choice == 1);
        y = results.EV(includeIdx & results.Choice == 1);
        scatter(x, y,'MarkerEdgeAlpha',0,'MarkerFaceColor',[255 42 116]/255); % risky choice

        x = results.ExpTrial(includeIdx & results.Choice == 2);
        y = results.EV(includeIdx & results.Choice == 2);
        scatter(x, y,'MarkerEdgeAlpha',0,'MarkerFaceColor',[0 240 138]/255); % safe choice

        x = results.ExpTrial(~includeIdx);
        y = results.EV(~includeIdx);
        scatter(x, y, 75, 'x','MarkerEdgeColor',[.5 .5 .5],'LineWidth',2);
        
        % Format
        plot(ax.XLim,[1 1],'Color',[0 0 0]);
        set(gca,'TickLength',[0 0])
        xlabel('Trials')
        ylabel('EV')
        title([num2str(subject) ': Accuracy = ' num2str(round(nanmean(results.Acc)*100,2)) '%'])
        
        set(gcf,'Position',[124 643 1667 468])
        
        % %%%%%%%
        % %%%%%%% Logistic function between EV and choice
        H2 = figure;

        x1 = results.EV(includeIdx & results.Choice == 1);
        y1 = ones(length(x1),1);
        scatter(x1,y1,'MarkerEdgeAlpha',0,'MarkerFaceColor',[255 42 116]/255); hold on % risky

        x2 = results.EV(includeIdx & results.Choice == 2);
        y2 = zeros(length(x2),1);
        scatter(x2,y2,'MarkerEdgeAlpha',0,'MarkerFaceColor',[0 240 138]/255); hold on % safe

        x = unique(results.EV(includeIdx));
        y1 = nan(length(x),1);
        y2 = y1;
        for i = 1:length(x)
            y1(i) = sum(results.Choice(includeIdx & results.EV == x(i)) == 1);
            y2(i) = sum(results.Choice(includeIdx & results.EV == x(i)) == 2);
        end

        [x, sortidx] = sort(results.EV(includeIdx));
        y = results.Choice(includeIdx)-1; % binary (0, 1) for choices (risky, safe)
        y = y(sortidx);
        [B,dev,stats] = glmfit(x,y,'binomial','logit');
        mfit = glmval(B,x,'logit');
        plot(x,1-mfit,'Color',[0 146 240]/255,'LineWidth',2);

        ax = gca;
        plot(ax.XLim,[.5 .5],'--k'); % 50% chance line

        ylim([ax.YLim(1)-.1 ax.YLim(2)+.1])
        xlabel('EV')
        ylabel('Risk-taking probability')

        IP = mfit(find(abs(x - 1) == min(abs(x - 1)))); % indifference point
        IP = IP(1);
        plot([IP IP],ax.YLim,'Color',[0 146 240]/255,'LineStyle',':','LineWidth',1.5);
        set(gca,'TickLength',[0 0])

        title([num2str(subject) ': Indifference point = ' num2str(round(IP,2))])
        set(gcf,'Position',[368 435 798 620])
        
    end
    
end

save(fullfile(dir_results,'groupBehaviour.mat'),'data'); % wide format (.mat)

writetable(...
    array2table([str2double(subjects) data.questionnaires],'variablenames',...
    {'Subject','IUS','Worry','Ethical','Financial','Health','Recreational','Social'}),...
    fullfile(dir_results,'groupBehaviour.csv')); % wide format (CSV)

writetable(T,fullfile(dir_results,'longData.csv')); % long format
