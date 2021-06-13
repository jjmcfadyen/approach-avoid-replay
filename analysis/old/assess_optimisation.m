%% Check how well the lambdas were optimised for sequenceness

clear all
clc

nType = 1; % null data
e = 1; % 1 = planning, 2 = transition display

loadSeq = false;

%% Directories

addpath('utils')
[sinfo, dinfo] = dir_cfg();

addpath(dinfo.tb_fieldtrip);
ft_defaults;

dir_load = fullfile(dinfo.data_meg,'seqOptimisation',['1vsRest_n' num2str(nType)]);
dir_save = fullfile(dinfo.results_meg,'sequenceness');

%% Parameters

subjects = sinfo.subjects;
N = length(subjects);

lags = 10:10:600; % in ms

timebins = 100:10:300; % in ms
orig_timebins = 0:10:300; % in ms
nT = length(timebins);

epochs = {'decision','transition','image','outcome'};
nE = length(epochs);

nLambda = 100;

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
U = U(uIdx,:);

nPerm = size(U,1);

if e == 1
    removeForced = true;
elseif e == 2
    removeForced = false;
end

%% Set parameters for optimisation

sidx = setdiff(1:length(subjects),[]);
N = length(sidx);

byParticipant = true; % maximise PARTICIPANT wave or GROUP wave
averageTime   = false;

by        = 'either'; % 'either' (fwd OR bwd), 'diff' (fwd-bwd), 'any' (fwd OR bwd OR fwd-bwd)
showPlots = false;

opWindow     = [0 300]; % restrict window within which to look
opWindow_idx = lags >= opWindow(1) & lags <= opWindow(2);

subtractNull = false;

%% Load data

load(fullfile(dinfo.results_meg_classifiers,['1vsRest_n' num2str(nType)],'CA_subject.mat')); % loads 'Y_subj' variable
CA = squeeze(mean(Y_subj(:,ismember([0:10:300],timebins),:),3)); % classifier accuracy

load(fullfile(dinfo.results_meg_classifiers,['1vsRest_n' num2str(nType)],'bestLambdas.mat')); % loads 'bL' variable
cBL = squeeze(bL(:,ismember([0:10:300],timebins),end)); % best lambdas per classifier training time

bestTime = find(mean(CA) == max(mean(CA)));
timerange = find(mean(CA) > quantile(mean(CA),.8));

nTW = length(timerange);
stimewin = nan(N,1);

%% Generate sequenceness

% Start with lambdas set to best classifier accuracy
SQ = nan(N,nTW,nPerm,3,length(lags));
grand = nan(N,nPerm,3,length(lags));
for s = 1:N

    subject = subjects{s};
    
    % Get epoch data to build sequenceness
    [fl,event,in] = pp_cfg(subject,'Task');
    [test, ~] = getTask(subject,fl,event,epochs{e});
    in.Fs = test.Fs;
    nTrls = size(test.L,1);

    if removeForced
        idx = event.table.Forced(ismember(event.table.ExpTrial,test.L(:,1))) == 0;
        test.D = test.D(idx);
        test.L = test.L(idx,:);
        test.x = test.x(idx);
    end
    nTrls = length(test.D);

    for tw = 1:nTW
        
        % Get classification accuracy & best lambda
        this_tw = timerange(tw);
        load(fullfile(dinfo.data_meg_classifiers,['1vsRest_n' num2str(nType)],subject,...
            [subject '_cc_train-' num2str(timebins(this_tw)) 'ms.mat'])); % loads 'output' variable
        this_lambda = cBL(s,this_tw);
        
        % Build sequenceness using this classifier
        classifier = cc_build(subject,timebins(this_tw),1,'1vsRest');
        classifier.B = classifier.B(:,:,this_lambda);
        classifier.I = classifier.I(:,this_lambda)';
        
        sq = [];
        for trl = 1:nTrls
            pred = cc_predict(test.D{trl},classifier,1,1); 
            sq(trl,:,:,:,:) = ss_build(pred,U,in);
        end
        sq = squeeze(mean(mean(sq,1),4));
        
        % insert into variable
        SQ(s,tw,:,:,:) = sq;
    end
end

grand = squeeze(SQ(:,timerange==bestTime,:,:,:));

% Plot
opts = [];
opts.x = lags;
opts.showOverall = true;
opts.showError = true;
opts.subtractNull = false;
opts.avCol = [0 0 1];
figure
for g = 1:3
    subplot(1,3,g)
    opts.g = g;
    ss_plot(squeeze(grand(:,:,g,:)),opts);
end
orig_grand = grand;

%% Find which training time gives best sequenceness signal per participant

if ~averageTime

    % Alter each subject's lambda and assess how it changes the group result
    for s = 1:N

        subject = subjects{sidx(s)};

        criteria = [];
        if byParticipant
            allTraces = nan(nTW,nPerm,3,length(lags));
        else
            allTraces = nan(nTW,N,nPerm,3,length(lags));
        end
        for tw = 1:nTW

            l = cBL(s,timerange(tw));

            if byParticipant
                tmp = squeeze(SQ(sidx(s),tw,:,:,:));
                allTraces(tw,:,:,:) = tmp;
            else
                tmp = grand;
                tmp(s,:,:,:) = squeeze(SQ(sidx(s),tw,:,:,:));
                allTraces(tw,:,:,:,:) = tmp;
            end

            criteria.timewindow(tw,1) = tw;
            criteria.lambda(tw,1) = l;
        end

        % assess peak changes per permutation  
        if byParticipant

           if subtractNull
                pGrand = [];
                pTrace = [];
                for g = 1:3
                    
                    thisg = squeeze(allTraces(timerange==bestTime,:,g,opWindow_idx));
                    thisp = squeeze(allTraces(:,:,g,opWindow_idx));
                    if g == 3
                        pGrand(1,g,:) = abs(thisg(1,:)) - quantile(max(abs(thisg(2:end,:)),[],2),.975);
                        for tw = 1:nTW
                            np = squeeze(thisp(tw,2:end,:));
                            pTrace(tw,g,:) = abs(squeeze(thisp(tw,1,:))) - quantile(max(abs(np),[],2),.975);
                        end
                    else
                        pGrand(1,g,:) = thisg(1,:) - quantile(max(thisg(2:end),[],2),.95);
                        for tw = 1:nTW
                            np = squeeze(thisp(tw,2:end,:));
                            pTrace(tw,g,:) = squeeze(thisp(tw,1,:)) - quantile(max(np,[],2),.95);
                        end
                    end
                end
            else
                pGrand = [];
                pGrand(1,:,:) = squeeze(allTraces(timerange==bestTime,1,:,opWindow_idx));
                pTrace = squeeze(allTraces(:,1,:,opWindow_idx));
            end

            [M,I] = max(pTrace,[],3);

            tmp = [];
            tmp.peakLoc = I;
            tmp.peakVal = M;

            for g = 1:3
                tmp.meanVal(:,g) = max(squeeze(pTrace(:,g,:)),[],2) - min(squeeze(pTrace(:,g,:)),[],2);
            end

            tmp.startVal = squeeze(pTrace(:,:,1));

            tmp.startVal(:,3) = abs(tmp.startVal(:,3));
            tmp.peakVal(:,3) = abs(tmp.peakVal(:,3));

        else
            pTrace = squeeze(allTraces(:,:,1,:,opWindow_idx));
            pGrand = squeeze(grand(:,1,:,opWindow_idx));
            tmp = so_criteria(pTrace,pGrand,opts);  
        end

        switch by
            case 'either'
                peakLoc = squeeze(round(max(tmp.peakLoc(:,1:2),[],2)));
            case 'diff'
                peakLoc = squeeze(round(tmp.peakLoc(:,3)));
            case 'any'
                peakLoc = squeeze(round(max(tmp.peakLoc(:,1:3),[],2)));
        end
        meanVal = tmp.meanVal;
        peakVal = tmp.peakVal;
        startVal = tmp.startVal;

        oidx = [1:size(peakLoc,1)]'; 

        % get best lambda
        switch by
            case 'either'
                cv = max(peakVal(:,1:2),[],2); % max of fwd OR bwd
            case 'any'
                cv = max(peakVal(:,1:3),[],2); % max of fwd OR bwd OR fwd-bwd
            case 'diff'
                cv = peakVal(:,3); % max fwd-bwd
        end
        tmp = sortrows([oidx cv(oidx) startVal(oidx) < cv(oidx)],[2 3],'descend');

        I = tmp(1,1);
        mv = tmp(1,2);

        mtw = criteria.timewindow(I,1);
        mp  = criteria.lambda(I,1);

        performance = [mv mp mtw]; % peak value change, lambda, time

        %{
        % plot lambdas for this permutation
        cmap = colours(size(pTrace,1),'plasma');
        figure
        for g = 1:3
            subplot(1,3,g)
            for tw = 1:nTW
                if byParticipant
                    y = squeeze(mean(pTrace(tw,g,:),2));
                else
                    y = squeeze(mean(pTrace(tw,:,g,:),2));
                end
                if outliers(tw)
                    plot(lags(opWindow_idx),y,'color',cmap(tw,:),'linewidth',1.2,'linestyle',':'); hold on
                elseif this_twin(tw) == mtw
                    plot(lags(opWindow_idx),y,'color',cmap(tw,:),'linewidth',2); hold on
                else
                    plot(lags(opWindow_idx),y,'color',cmap(tw,:)); hold on
                end
            end
            set(gca,'ticklength',[0 0])
            ax = gca;
            plot(ax.XLim,[0 0],'k:'); hold on
        end
        %}

        % pick best lambda per permutation
        tmp = grand;
        tmp(s,:,:,:) = squeeze(SQ(sidx(s),... % subject
                               performance(3),... % time
                               :,... % permutation
                               :,:)); % g, lags

        if showPlots

            opts = [];
            opts.x = lags;
            opts.showOverall = true;
            opts.showError = true;
            opts.subtractNull = true;

            figure;
            cmap = [0 0 1; 1 0 0];
            for i = 1:2
                opts.avCol = cmap(i,:);
                if i == 1
                    y = grand;
                elseif i == 2
                    y = tmp;
                end
                for g = 1:3
                    subplot(1,3,g)
                    opts.g = g;
                    ss_plot(squeeze(y(:,:,g,:)),opts);
                end
            end
            sgtitle(subject)
            drawnow


            opts = [];
            opts.x = lags;
            opts.showOverall = true;
            opts.showError = true;
            opts.nullPerms_overall = true;
            opts.nullMoving = false;

            gt = squeeze(mean(tmp));
            gg = squeeze(mean(grand));
            ylims = [min([min(gt(:)) min(gg(:))]) max([max(gt(:)) max(gg(:))])];

            for i = 1:2
                figure
                if i == 1
                    thisd = grand;
                    ttext = 'Before';
                elseif i == 2
                    thisd = tmp;
                    ttext = 'After';
                end
                for g = 1:3
                    subplot(1,3,g)
                    opts.g = g;
                    ss_plot(squeeze(thisd(:,:,g,:)),opts);
                    ylim(ylims)
                end
                sgtitle(ttext);
            end

        end

        % save to variable & update grand
        stimemax(s,1) = performance(1);
        stimewin(s,1) = performance(3);
        bL(s,bestTime,:) = performance(2);
        grand(s,:,:,:) = tmp(s,:,:,:);
    end

else

    for s = 1:N
        tmp = [];
        for tw = 1:nTW
            tmp(tw,:,:,:) = squeeze(SQ(s,timerange(tw),CA(s,timerange(tw),2),:,:,:));
        end
        grand(s,:,:,:) = squeeze(mean(tmp));
    end

end

% plot final    
opts = [];
opts.x = lags;
opts.showOverall = true;
opts.showError = true;
opts.subtractNull = false;
opts.avCol = [1 0 0];
%     figure
for g = 1:3
    subplot(1,3,g)
    opts.g = g;
    ss_plot(squeeze(grand(:,:,g,:)),opts);
end
sgtitle([num2str(timebins(bestTime)) ' ms'])
drawnow


save(fullfile(dir_save,['seq_bestLambdas_' epochs{e} '.mat']),'cBL','CA','stimewin','stimemax');
