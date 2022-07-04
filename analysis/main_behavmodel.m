%% Behavioural modelling

clear all
clc

%% Directories & parameters

cd D:\2020_RiskyReplay\approach-avoid-replay
addpath('utils');

dir_data = 'D:\2020_RiskyReplay\data\behav';
dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_output = 'D:\2020_RiskyReplay\results\modelling';

parameters = get_parameters(dir_raw);

% get subjects
subjects = unique(parameters.schar);
N = length(subjects);

%% Simulate models to get performance

models = set_models();
nModels = length(models);
mnames = extractfield(models,'name');

rng(1);

numStart = 1;

% get unique datasets
tmp = [];
sidx = [];
for s = 1:N
    disp(['Reading in subject ' num2str(s) ' data...'])
    d = parse_behav(subjects{s},dir_data); 
    try
        tmp = [tmp; d.EV'];
        sidx(s) = s;
    end
end
[~,~,uD] = unique(tmp,'rows');

nIterations = length(unique(uD));

datasets = [];
for it = 1:nIterations
    datasets(it) = find(uD==it,1,'first');
end

optfit = [];
optfit.nLL = nan(nIterations,nModels);
optfit.acc = nan(nIterations,nModels);
optfit.params = nan(nIterations,nModels,10);
for it = 1:2%nIterations

    disp('====================================================================================')
    disp(['ITERATION ' num2str(it)])
    disp('====================================================================================')

    d = parse_behav(subjects{datasets(it)},dir_data); % randomly select subject protocol for simulations
    
    d.Choice = nan(size(d,1),1); % make choices perfect
    d.Choice(d.EV>1 & d.Forced==0) = 1;
    d.Choice(d.EV<1 & d.Forced==0) = 0;
    d.Choice(d.EV==1 & d.Forced==0) = rand(sum(d.EV==1 & d.Forced==0))>.5;

    d.Transition(d.Choice==0) = 0; % update transitions
    d.Outcome(d.Choice==0) = 1;
    idx = find(d.Choice==1 & d.Forced==0);
    for i = 1:length(idx)
        d.Transition(idx(i)) = datasample([1 2],1,'replace',false,'weights',[d.P(idx(i)) 1-d.P(idx(i))]);
        d.Outcome(idx(i)) = d.(['nV_' num2str(d.Transition(idx(i)))])(idx(i));
    end
    
    % Get performance of each model
    nLL = nan(nModels,1);
    BIC = nan(nModels,1);
    ACC = nan(nModels,1);
    params = nan(nModels,size(optfit.params,3));
    parfor m = 1:nModels
        
        disp(['Model ' num2str(m) '...'])

        % calculate accuracy of predictions
        idx = d.Forced==0 & d.EV~=1;
        if contains(mnames{m},'random')
            disp(['Running random iterations for dataset ' num2str(it) ', model ' num2str(m) '...'])
            acc = nan(1,100);
            err = nan(1,100);
            bic = nan(1,100);
            thisparams = [];
            for i = 1:100
                warning off
                [~,fitvals,~] = fit_startrange(models(m),d,numStart,'one');
                warning on;
                [err(i), output] = sim_model(d,models(m),fitvals','one');
                if isnan(err(i))
                    error('Nan in nLL')
                end
                acc(i) = mean(output.T.y(idx) == output.T.yhat(idx));
                [~,bic(i)] = aicbic(err(i) * (-1),length(fitvals),sum(d.Forced==0));
                thisparams(i,:) = fitvals;
            end
            err = mean(err);
            bic = mean(bic);
            params(m,:) = [mean(thisparams) nan(1,size(optfit.params,3)-size(thisparams,2))];
        else
            warning off
            [startvals,fitvals,err] = fit_startrange(models(m),d,numStart,'both');
            warning on

            beststart = find(err==min(err),1,'first');
            x = fitvals(:,beststart)';

            [err, output] = sim_model(d,models(m),x,'both');
            [aic,bic] = aicbic(err * (-1),length(x),sum(d.Forced==0));
            acc = output.T.y(idx) == output.T.yhat(idx);

            if length(fitvals) < size(optfit.params,3)
                x = [x nan(1,size(optfit.params,3)-length(x))];
            end
            params(m,:) = x;
        end

        if any(isnan(acc)) || any(isnan(bic)) || any(isnan(err)) || any(isinf(acc)) || any(isinf(bic)) || any(isinf(err))
            error('NaNs in accuracy/bic/nLL');
        end

        % log results
        nLL(m,:) = err;
        BIC(m,:) = bic;
        ACC(m,:) = mean(acc);

        disp(['Iteration ' num2str(it) ', model: ' models(m).name ' ------ accuracy = ' num2str(round(mean(acc),2))]);

    end

    optfit.nLL(it,:) = nLL;
    optfit.bic(it,:) = BIC;
    optfit.acc(it,:) = ACC;
    optfit.params(it,:,:) = params;

end

% Save
save(fullfile(dir_output,'comparestrategies.mat'),'optfit','models');

% Plot
mnames = extractfield(models,'name');

modelOrder = [1 2 4 6 3 5 7:11 13 15 12 14 16:20 22 24 21 23 25:32];

figure
for i = 1:3

    subplot(1,3,i)

    if i==1
        y = optfit.nLL(1:2,modelOrder);
        thistitle = 'Negative Log Likelihood';
    elseif i==2
        y = optfit.bic(1:2,modelOrder);
        thistitle = 'BIC';
%         y = y(:,end) - y;
    elseif i==3
        y = optfit.acc(1:2,modelOrder);
        thistitle = 'Accuracy';
    end

    if size(y,1)>1
        if i==3
            m = mean(y);
        else
            m = sum(y);
        end
    else
        m=y;
    end
    if i==2
        m = max(m)-m;
    end
    err = std(y);%/sqrt(size(y,1));
    upper = m+err;
    lower = m-err;

    bar(m); hold on
    for j = 1:size(y,2)
        plot([j j],[upper(j) lower(j)],'k'); 
    end
    
    set(gca,'ticklength',[0 0])
    title(thistitle)
    set(gca,'xtick',1:nModels)
    set(gca,'xticklabels',mnames(modelOrder));

end

%% Fit models

models = set_models();
nModels = length(models);

onCluster = true;
makePlots = false;
numStart = 1; % number of different starting values to use per parameter

% create empty structures for model optimisation output
if ~onCluster
    optim = struct();
    for m = 1:nModels

        np = 0;
        freeidx = [];
        cc = 0;
        for p = fieldnames(models(m).params)'
            cc = cc+1;
            if isnan(models(m).params.(p{1}).val)
                np = np+1;
                freeidx = [freeidx cc];
            end
        end
        pnames = fieldnames(models(m).params);
        pnames = pnames(freeidx);

        optim(m).nLL = nan(N,1); % negative log likelihood
        optim(m).params = array2table(nan(N,np),'variablenames',pnames');
        optim(m).aic = nan(1,N);
        optim(m).bic = nan(1,N);
        optim(m).acc = nan(1,N);
        optim(m).startvals.start = nan(N,np,numStart^np); % subjects, parameters, starting combination index
        optim(m).startvals.fit = nan(N,np,numStart^np); % subjects, parameters, starting combination index
        optim(m).startvals.nLL = nan(N,numStart^np);
        optim(m).model = models(m);
    end
end

% loop through subjects
for s = 1:N
   
    % Read in data
    d = parse_behav(subjects{s},dir_data);
    
    % Exclude missed trials
    d = d(d.RT<30,:);

    % Convert choice 2 (avoid) to 0
    d.Choice(d.Choice==2) = 0;

    % Ignore forced-choice trials
    idx = d.Forced==0;

    if makePlots
        plotdata = [];
    end

    if onCluster
        dir_save = fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',subjects{s});
        if ~exist(dir_save)
            mkdir(dir_save)
        end
        cc = 0;
    end

    for m = 1:nModels

        if ~onCluster

            [startvals,fitvals,err] = fit_startrange(models(m),d,numStart,'both');
    
            optim(m).startvals.start(s,:,:) = startvals;
            optim(m).startvals.nLL(s,:) = err;
            optim(m).startvals.fit(s,:,:) = fitvals;
    
            beststart = find(err==min(err));
            x = mean(fitvals(:,beststart),2);
    
            % predict data using optimised parameters
            [~,output] = sim_model(d,models(m),x,'both');
    
            if makePlots
                
                mdl_y = fitglm(output.T(idx,:),'y~EV','distribution','binomial');
                mdl_yhat = fitglm(output.T(idx,:),'yhat~EV','distribution','binomial');
     
                pred_y = predict(mdl_y);
                pred_yhat = predict(mdl_yhat);
    
                [sorted,sortidx] = sort(output.T.EV(idx));
                
                y = output.T.y(idx);
                yhat = output.T.yhat(idx);
    
                plotdata.x{m} = sorted;
                plotdata.y{m} = y(sortidx);
                plotdata.yfit{m} = pred_y(sortidx);
                plotdata.yhat{m} = yhat(sortidx);
                plotdata.yhatfit{m} = pred_yhat(sortidx);
            end
    
            [aic,bic] = aicbic(-output.nLL,length(x),size(output.T,1));
            
            optim(m).nLL(s,:) = output.nLL;
            optim(m).params(s,:) = array2table(x');
            optim(m).aic(s) = aic;
            optim(m).bic(s) = bic;
            optim(m).acc(s) = mean(output.T.y(idx)==output.T.yhat(idx));
        else
            
            cc = cc + 1;
            filename = fullfile(dir_save,['modellingbatch-it' num2str(cc) '.mat']);
            
            model = models(m);
            if isempty(model)
                error('Empty model')
            end

            info = [];
            info.subject = subjects{s};
            info.m = m;
            info.randtype = 'both';

            thisd = d;

            save(filename,'model','thisd','numStart','info');
        end
    end

    if onCluster
        filename = fullfile(dir_save,['modellingbatch-it$SGE_TASK_ID.mat']);
        generate_jobs_modelling(filename,'holly',cc);
    end

    % plot logistic regression slopes for actual data vs model-predicted data
    if makePlots && ~onCluster
    
        figure
        sgtitle(['Subject ' num2str(s) ': ' subjects{s}]);

        for m = 1:nModels

            subplot(3,4,m)
            scatter(plotdata.x{m},plotdata.y{m},'markerfacecolor','k','markeredgecolor','none','markerfacealpha',.5); hold on
            scatter(plotdata.x{m},plotdata.yhat{m},'r'); hold on
            
            plot(plotdata.x{m},plotdata.yfit{m},'k'); hold on
            plot(plotdata.x{m},plotdata.yhatfit{m},'r'); hold on
    
            set(gca,'ticklength',[0 0])
            title([models(m).name ': ' num2str(round(mean(plotdata.y{m}==plotdata.yhat{m})*100,2)) '%, BIC=' num2str(round(optim(m).bic(s)))])
            xlabel('EV')
            ylabel('Choice probability')
        end
        set(gcf,'position',[1 31 1920 973])
        drawnow
    end
end

% Save output
if ~onCluster
    save(fullfile(dir_output,'modelfitting_randtype-both.mat'),'models','optim');
end

%% Plot

load(fullfile(dir_output,'modelfitting.mat'));

excludedSubjects = {}; %{'263098','680913'};
includedidx = ~ismember(subjects,excludedSubjects);

% get model names
mnames = cell(1,length(optim));
for m = 1:length(mnames)
    mnames{m} = optim(m).model.name;
end

% get order in model order for plotting
modelOrder = [1 12 27 29 23 28 30 31 32 2 3 5 7 4 6 8 9 10 11 13 15 17 14 16 18 19 20 21 22 24 25 26];

% get negative log likelihoods
nLL = nan(N,nModels);
for m = 1:nModels
    nLL(:,m) = optim(m).nLL;
end

% calculate AIC and BIC
AIC = nan(N,nModels);
BIC = nan(N,nModels);
ACC = nan(N,nModels);
for m = 1:nModels
    AIC(:,m) = optim(m).aic;
    BIC(:,m) = optim(m).bic;
    ACC(:,m) = optim(m).acc;
end

% group model evidence
[sorted,sortidx] = sort(sum(BIC));

figure
bar(sorted)
set(gca,'xtick',1:nModels); set(gca,'xticklabels',mnames(sortidx));
title('group model evidence')

% subject winning model histogram
[~,winner] = min(BIC,[],2);

[sorted,sortidx] = sort(BIC,2);
windiff = sorted(:,2)-sorted(:,1);
for s = 1:N
%     if windiff(s) < 3
%         toprank = sortidx(s,sorted(s,:) < sorted(s,1)+3);
%         if any(toprank==1)
%             winner(s) = 1;
%         end
%         if winner(s) ~=1 && any(contains(mnames(toprank),'twopath'))
%             winner(s) = toprank(contains(mnames(toprank),'twopath'));
%         end
%     end
    winner(s) = find(modelOrder==winner(s));
end

orderednames = mnames(modelOrder);

figure
histogram(winner,nModels)
set(gca,'xtick',unique(winner)); set(gca,'xticklabels',orderednames(unique(winner)));
title('Winning model per subject')

tmp = [];
for m = 1:nModels
    tmp(m,:) = [m sum(winner==m) ];
end
tmp = array2table(tmp(tmp(:,2)>0,:),'variablenames',{'model','count'});

tmp.calculation = cell(size(tmp,1),1);
tmp.pathtype = cell(size(tmp,1),1);
tmp.qtype = cell(size(tmp,1),1);
for m = 1:size(tmp,1)
    tmp.calculation{m} = optim(strcmp(mnames,orderednames{tmp.model(m)})).model.calctype;
    tmp.pathtype{m} = optim(strcmp(mnames,orderednames{tmp.model(m)})).model.pathtype;
    tmp.qtype{m} = optim(strcmp(mnames,orderednames{tmp.model(m)})).model.qtype;
end

% make table to put into replay analyses
winnertable = array2table(subjects,'variablenames',{'Subject'});
winnertable.bModelNum = winner;

winnertable.bCalculationType = cell(N,1);
winnertable.bPathChoice = cell(N,1);
for s = 1:N
    winnertable.bCalculationType{s} = tmp.calculation{tmp.model==winner(s)};
    winnertable.bPathChoice{s} = tmp.pathtype{tmp.model==winner(s)};
end

writetable(winnertable,'D:\2020_RiskyReplay\results\modelling\behavmodeltable.csv');

% get acc
acc = nan(N,1);
for s = 1:N

    d = parse_behav(subjects{s},dir_data);
    d = d(d.RT<30 & d.Forced==0 & d.EV~=1,:);

    d.acc = (d.EV>1 & d.Choice==1) | (d.EV<1 & d.Choice==2);
    acc(s,1) = mean(d.acc);

end

% compare accuracy between caching and calculating
gacc = cell(1,2);
for g = 1:2
    if g==1
        sidx = find(contains(winnertable.bCalculationType,'calculating'));
    else
        sidx = find(contains(winnertable.bCalculationType,'caching'));
    end
    gacc{g} = acc(sidx);
end

[h,p,~,stats] = ttest2(gacc{1},gacc{2});

% % show individual subjects, ordered by accuracy
% [sorted,sortidx] = sort(acc,'descend');
% 
% figure
% for s = 1:N
%     subplot(4,7,s)
%     bic = BIC(sortidx(s),:);
%     bic = bic(strcmp(mnames,'null')) - bic(:);
%     bic = bic(modelorder);
% %     tmp = [[2:10]' bic(2:end)];
%     bar(bic);
%     title([num2str(round(acc(sortidx(s))*100)) '% - ' mnames{find(bic==max(bic))}])
%     set(gca,'ticklength',[0 0])
% end

% parameter estimates
figure
for m = 1:nModels

    subplot(8,8,m)
    
    paramnames = fieldnames(models(m).params);
    np = size(optim(m).params,2);
    for p = 1:np
        
        y = table2array(optim(m).params(includedidx,p));

        scatter(ones(sum(includedidx),1)*p,y); hold on
        scatter(p,median(y),'markerfacecolor','k'); hold on

    end   
    xlim([0 np+1])
    set(gca,'xtick',1:np)
    set(gca,'xticklabels',paramnames)
    title(models(m).name)

end

%% Parameter recovery & Model generalisability

models = set_models();
nModel = length(models);

onCluster = true;

nIterations = 50;
numStart = 1;

maxnp = 0;
for m1 = 1:nModels
    pnames = fieldnames(models(m1).params);
    np = length(pnames);
    freeidx = zeros(1,np);
    nFree = 0;
    for p = 1:np
        if isnan(models(m1).params.(pnames{p}).val)
            freeidx(p) = 1;
            nFree = nFree + 1;
        end
    end
    if nFree > maxnp
        maxnp = nFree;
    end
end

% run procedure
if ~onCluster
    nLL = nan(nModels,nModels,nIterations);
    AIC = nan(nModels,nModels,nIterations);
    BIC = nan(nModels,nModels,nIterations);
    recovery = cell(1,nModels);
    paramlog = nan(nModels,nModels,nIterations,maxnp);
end

% get unique datasets
tmp = [];
sidx = [];
for s = 1:N
    disp(['Reading in subject ' num2str(s) ' data...'])
    d = parse_behav(subjects{s},dir_data); 
    try
        tmp = [tmp; d.EV'];
        sidx(s) = s;
    end
end
[~,~,uD] = unique(tmp,'rows');

datasets = [];
for it = 1:length(unique(uD))
    datasets(it) = find(uD==it,1,'first');
end

D = cell(1,2);
for it = 1:2 % last iteration only applied to 2 subjects
    D{it} = parse_behav(subjects{find(uD==it,1,'first')},dir_data);
    D{it} = table2struct(D{it});
    D{it} = rmfield(D{it},'Choice');
    D{it} = rmfield(D{it},'Transition');
    D{it} = rmfield(D{it},'Outcome');
    D{it} = rmfield(D{it},'RT');
    D{it} = rmfield(D{it},'Timestamp');
    D{it} = rmfield(D{it},'Acc');
    D{it} = rmfield(D{it},'bAcc');
    D{it} = rmfield(D{it},'Subject');
    D{it} = struct2table(D{it});
end

for m1 = 1:nModels

    disp('=====================================================')
    disp(['MODEL ' num2str(m1)])
    disp('=====================================================')

    pnames = fieldnames(models(m1).params);
    np = length(pnames);
    freeidx = zeros(1,np);
    nFree = 0;
    for p = 1:np
        if isnan(models(m1).params.(pnames{p}).val)
            freeidx(p) = 1;
            nFree = nFree + 1;
        end
    end
    freeidx = find(freeidx);

    if onCluster
        batchdir = fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',['m' num2str(m1)]);
        if ~exist(batchdir)
            mkdir(batchdir);
        end
        cc = 0;
    end

    % Find corresponding model in 'optim', if using subject parameters to constrain random generation
    optimModelNum = NaN;
    for i = 1:nModels
        if strcmp(optim(i).model.name,models(m1).name)
            optimModelNum = i;
        end
    end

    % Generate choices using lots of different sets of parameters
    origP = nan(nIterations,nFree); % 1 = iteration number, 2 = parameter number
    for it = 1:nIterations
     
        % generate parameter set
        itparams = nan(nFree,1);
        for p = 1:nFree
            
%             fitparams = table2array(optim(optimModelNum).params(:,p));
%             thisparam = normrnd(mean(fitparams),std(fitparams));
%             thisrange = [min(fitparams) max(fitparams)];
%             if isinf(thisrange(2))
%                 thisrange(2) = 10; % set hard limit - prevents NaNs created from high tau, for example
%             end
%             while thisparam < thisrange(1) || thisparam > thisrange(2)
%                 thisparam = normrnd(mean(fitparams),std(fitparams)); % make sure the randomly generated number is within the correct range
%             end
% %             thisparam = unifrnd(min(fitparams),max(fitparams),1)

            switch pnames{p}
                case 'threshold'
                    thispdf = truncate(makedist('normal',0,6),-12,12);
                    thisparam = random(thispdf,1,1);
                case 'rewthreshold'
                    thispdf = truncate(makedist('normal',0,6),0,12);
                    thisparam = random(thispdf,1,1);
                case 'lossthreshold'
                    thispdf = truncate(makedist('normal',0,6),-12,0);
                    thisparam = random(thispdf,1,1);
                case 'tau'
                    thisparam = datasample([0.1:0.1:5 5.5:0.5:10],1);
                case 'alpha'
                    thisparam = datasample(0:0.01:1,1);
                otherwise
                    error('Missing parameter name')
            end

            itparams(p,1) = thisparam;
        end
        origP(it,:) = itparams;
    end

    if ~onCluster

        newP = nan(nIterations,nFree); % 1 = iteration number, 2 = parameter number
        t0 = tic;
        parfor it = 1:nIterations
    
            % simulate new data
            d = D{randi(length(D))} % randomly select subject protocol for simulations

            % generate predicted data
            [~, output, thisd] = sim_model(d,models(m1),squeeze(origP(it,:)),'one');
    
            % predict data using each model
            for m2 = 1:nModels

                disp(['--------- M1 = ' num2str(m1) ', iteration = ' num2str(it) ', M2 = ' num2str(m2) '...'])

                [startvals,fitvals,err] = fit_startrange(models(m2),thisd,numStart,'both');
        
                beststart = find(err==min(err));
                x = mean(fitvals(:,beststart),2);
    
                if m2==m1
                    newP(it,:) = x;
                end
    
                [~,output] = sim_model(thisd,models(m2),x);
    
                [aic,bic] = aicbic(output.nLL*(-1), nFree, size(output.T,1));
                AIC(m1,m2,it) = aic;
                BIC(m1,m2,it) = bic;
                nLL(m1,m2,it) = output.nLL;
                paramlog(m1,m2,it,:) = [x' nan(1,maxnp-length(x))];
            end
        end
        timeleft = toc(t0); % in seconds

        if m1 < nModels
            disp(['Estimated time remaining = ' (nModels-m1)*timeleft/60 ' minutes.'])
        end
    else

        for it = 1:nIterations

            % simulate new data
            d = D{randi(length(D))};% randomly select subject protocol for simulations

            % generate predicted data
            [~, output, thisd] = sim_model(d,models(m1),squeeze(origP(it,:)),'one');

            % predict data using each model
            for m2 = 1:nModels

                cc = cc + 1;
                filename = fullfile(batchdir,['modellingbatch-it' num2str(cc) '.mat']);
    
                model = models(m2);
    
                info = [];
                info.m1 = m1;
                info.m2 = m2;
                info.it = it;
                info.origParams = squeeze(origP(it,:));
                info.nModels = nModels;
                info.nIterations = nIterations;
                info.maxnp = maxnp;
                info.randtype = 'both';
    
                save(filename,'model','numStart','thisd','info');
            end
        end
    end

    if ~onCluster
        P = origP;
        P(:,:,2) = newP;
    
        recovery{m1} = P;
    
        disp(['Model ' mnames{m1}])
        [~,idx] = min(squeeze(BIC(m1,:,:)));
        disp(['--- mean winning model = ' num2str(mean(idx==m1))])
        disp(['--- parameter recovery:'])
        corr(squeeze(P(:,:,1)),squeeze(P(:,:,2)))
    else
        filename = fullfile(batchdir,['modellingbatch-it$SGE_TASK_ID.mat']);
        generate_jobs_modelling(filename,'holly',cc);
    end
end

if ~onCluster
    save(fullfile(dir_output,'modelrecovery.mat'),'AIC','BIC','nLL','paramlog','d');
else
    warning('Run "cluster_organiseModels_recovery.m" on cluster to get everything in the right data structures')
end

%% Compute metrics on recovery

tmp = load(fullfile(dir_output,'modelrecovery_evendistribution.mat'));

modelorder = [1 2 4 6 3 5 7:11 13 15 12 14 16:20 22 24 21 23 25:32];

BIC = tmp.BIC(modelorder,modelorder,:);
recovery = tmp.recovery(modelorder);


% Plot parameter recovery
acc = nan(nModels,4);
sensitivity = nan(nModels,4);
specificity = nan(nModels,4);
for m1 = 1:nModels

    y = abs(corr(squeeze(recovery{modelorder(m1)}(:,:,1)),squeeze(recovery{modelorder(m1)}(:,:,2))));

    if any([1 2 3 4]==m1)
        pidx = [1 3];
    elseif any([5 6 7]==m1)
        pidx = [1 2 3];
    elseif any([8:13 17:22 26:27]==m1)
        pidx = [1 3 4];
    elseif any([14:16 23:25]==m1)
        pidx = [1:4];
    end

    for p = 1:size(y,2)

        TP = y(p,p);
        FP = mean(y(:,p));
        FN = mean(y(p,setdiff(1:size(y,2),p)));
        TN = mean(mean(y(setdiff(1:size(y,2),p),setdiff(1:size(y,2),p))));
    
        acc(m1,pidx(p)) = (TP+TN)/(TP+TN+FP+FN);
        sensitivity(m1,pidx(p)) = TP/(TP+FN);
        specificity(m1,pidx(p)) = TN/(TN+FP);
    end
end

figure
for m = 1:nModels

    thiscorr = corr(squeeze(recovery{modelorder(m)}(:,:,1)),squeeze(recovery{modelorder(m)}(:,:,2)));

    subplot(4,7,m)
    imagesc(abs(thiscorr))
    colormap('gray')
    caxis([0 1])
    title(models(modelorder(m)).name)
    set(gca,'xtick',1:size(thiscorr,2))
    set(gca,'xticklabels',models(modelorder(m)).paraminfo.names)
    set(gca,'ticklength',[0 0])

end

% Plot model specificity
y = nan(nModels,nModels);
for m1 = 1:nModels
    [~,idx] = min(squeeze(BIC(m1,:,:)));
    for m2 = 1:nModels
        y(m1,m2) = mean(idx==m2);
    end
end

acc = nan(1,nModels);
sensitivity = nan(1,nModels);
specificity = nan(1,nModels);
for m1 = 1:nModels
 
    TP = y(m1,m1);
    FP = sum(y(:,m1));
    FN = sum(y(m1,setdiff(1:nModels,m1)));
    TN = sum(sum(y(setdiff(1:nModels,m1),setdiff(1:nModels,m1))));

    acc(m1) = (TP+TN)/(TP+TN+FP+FN);
    sensitivity(m1) = TP/(TP+FN);
    specificity(m1) = TN/(TN+FP);

end

figure
imagesc(y')
caxis([0 1])
colormap('gray')
set(gca,'ticklength',[0 0])
set(gca,'xtick',1:length(modelorder))
set(gca,'xticklabels',mnames(modelorder))
set(gca,'ytick',1:length(modelorder))
set(gca,'yticklabels',mnames(modelorder))
xlabel('Generative model')
ylabel('Winning predictive model')

%% Optimisim/pessimism hybrid model - simulate

% D = parse_behav(subjects{11},dir_data); % use random subject's experiment protocol
% D = table2struct(D);
% D = rmfield(D,'Choice');
% D = rmfield(D,'Transition');
% D = rmfield(D,'Outcome');
% D = rmfield(D,'RT');
% D = struct2table(D);
% 
% modelcombos = [
%     1 55
%     1 19
%     1 37
%     19 19
%     37 37
%     ];
% 
% paramiterations = [];
% 
% paramiterations.alpha = [0.1 0.5 0.9];
% paramiterations.gain = [0.001 0.1 1 5];
% paramiterations.tau = [0.1 0.5 1 5];
% paramiterations.threshold = [1];
% paramiterations.forgetting = [0.1 0.5 0.9];
% 
% paramiterations.MB_forgetting = paramiterations.forgetting;
% paramiterations.MF_forgetting = paramiterations.forgetting;
% paramiterations.replayalpha = paramiterations.alpha;
% paramiterations.MB_threshold = paramiterations.threshold;
% 
% allparams = [];
% GAIN = [];
% hybridModels = [];
% for v = 1%:3 % change how volatile the probability changes are
% 
%     if v==1
%         d = D;
%     elseif v==2
%         d = changeVolatility(D,0.6);
%     elseif v==3
%         d = changeVolatility(D,0.3);
%     end
%     
%     cc = 0;
%     for m = 1:size(modelcombos,1)
%         for g = 1:2 % gain type (gain = MF-Hybrid, or Hybrid-MF)
% 
%             cc = cc+1;
%     
%             MB = modelcombos(m,1);
%             MF = modelcombos(m,2);
%         
%             % Merge models into hybrid model
%             hybrid = makeHybrid(models(MF),models(MB));
% 
%             if g==1
%                 hybrid.gaintype = 'gain';
%             elseif g==2
%                 hybrid.gaintype = 'gainneed';
%             end
% 
%             hybridModels{cc} = hybrid;
% 
%             pnames = fieldnames(hybrid.params);
%             np = length(pnames);
%             freeidx = [];
%             for p = 1:np
%                 if isnan(hybrid.params.(pnames{p}).val)
%                     freeidx = [freeidx p];
%                 end
%             end
%             pnames = pnames(freeidx);
%             np = length(pnames);
%         
%             % Randomly generate parameters
%             tmp = [];
%             for p = 1:np
%                 tmp{p} = paramiterations.(pnames{p});
%             end
%     
%             if np==4
%                 params = combvec(tmp{1},tmp{2},tmp{3},tmp{4})';
%             elseif np==5
%                 params = combvec(tmp{1},tmp{2},tmp{3},tmp{4},tmp{5})';
%             elseif np==6
%                 params = combvec(tmp{1},tmp{2},tmp{3},tmp{4},tmp{5},tmp{6})';
%             elseif np==7
%                 params = combvec(tmp{1},tmp{2},tmp{3},tmp{4},tmp{5},tmp{6},tmp{7})';
%             end
%             nIterations = size(params,1);
%         
%             allparams{cc} = params;
%         
%             % Generate data for each iteration
%             predreplay = nan(nIterations,4); % iteration, path/choice combination
%             gain = nan(nIterations,4); % iteration, path/choice combination
%             parfor it = 1:nIterations
%         
%                 disp(['MB = ' num2str(MB) ', MF = ' num2str(MF) ': Iteration ' num2str(it) ' of ' num2str(nIterations) '...'])
%                 [~, ~, thisd] = sim_model(d,hybrid,params(it,:),'one');
%         
%                 plotdata = [];
%                 plotdata(1) = mean(thisd.replayinitiated_rewarding(thisd.Choice==1));
%                 plotdata(2) = mean(thisd.replayinitiated_aversive(thisd.Choice==1));
%                 plotdata(3) = mean(thisd.replayinitiated_rewarding(thisd.Choice==0));
%                 plotdata(4) = mean(thisd.replayinitiated_aversive(thisd.Choice==0));
%                 predreplay(it,:) = plotdata;
% 
%                 plotdata = [];
%                 plotdata(1) = mean(thisd.gain_rewarding(thisd.Choice==1));
%                 plotdata(2) = mean(thisd.gain_aversive(thisd.Choice==1));
%                 plotdata(3) = mean(thisd.gain_rewarding(thisd.Choice==0));
%                 plotdata(4) = mean(thisd.gain_aversive(thisd.Choice==0));
%                 gain(it,:) = plotdata;
%         
%             end
%         
%             GAIN{v,cc}.replayinitiated = predreplay;
%             GAIN{v,cc}.gain = gain;
%         end
%     end
% end
% 
% % Plot
% v = 1;
% for m = 1:size(GAIN,2)
%   
%     y = GAIN{v,m}.gain;
%     nIterations = size(y,1);
% 
%     counterfactual = nan(nIterations,2);
%     counterfactual(:,1) = squeeze(y(:,2)-y(:,1)); % approach: aversive - rewarding
%     counterfactual(:,2) = squeeze(y(:,3)-y(:,4)); % avoid: rewarding - aversive
% 
%     hybrid = hybridModels{m};
%     pnames = fieldnames(hybrid.params);
%     np = size(allparams{m},2);
% 
%     figure;
%     for p = 1:np
%         subplot(2,ceil(np/2),p)
%         [sorted,sortidx] = sort(allparams{m}(:,p));
%         ii = 0;
%         for i = 1:2
%             plot(sorted,counterfactual(sortidx,i)); hold on
%         end
%         plot(sorted([1 end]),[0 0],'k:'); hold on
%         xlabel(pnames{p})
%         ylabel('Gain')
%     end
%     sgtitle(replace([hybrid.MB.name '   ' hybrid.MF.name],'_',' '))
% end


%% Optimisim/pessimism hybrid model - fit to subjects

error('Write something here - only select two-path strategy models? or choose best model and compare to all hybrids?')

onCluster = false;
numStart = 1;

replay = readtable('D:\2020_RiskyReplay\results\replay\replay_lme_planning.csv');

% Choose which model-based and model-free algorithms to merge into hybrid model
if ~onCluster
    BIC = nan(N,3);
    nLL = nan(N,3);
    ACC = nan(N,3);
end
allmodels = [];

modelcombos = [
    1 55   % how do mental calculations influence habitual actions?
    1 19   % how do mental calculations influence learning paths from experience?
    1 37   % how do mental calculations influence learning paths from experience (considering neg state)?
    19 19  % how does taking into account probability influence learning paths from experience?
    37 37  % how does taking into account probability influence learning paths from experience (considering neg state)?
    ];

allparams = [];
cc = 0;
for m = 1:size(modelcombos,1)

    MB = modelcombos(m,1);
    MF = modelcombos(m,2);

    hybrid = makeHybrid(models(MF),models(MB));

    pnames = fieldnames(hybrid.params);
    np = length(pnames);
    freeidx = [];
    for p = 1:np
        if isnan(hybrid.params.(pnames{p}).val)
            freeidx = [freeidx p];
        end
    end

    % Fix model-based threshold to 1
    tmp = fieldnames(hybrid.params);
    idx = find(contains(tmp,'threshold'));
    for p = 1:length(idx)
        hybrid.params.(tmp{idx(p)}).val = 1;
    end
    tmp = fieldnames(hybrid.MB.params);
    idx = find(contains(tmp,'threshold'));
    for p = 1:length(idx)
        hybrid.MB.params.(tmp{idx(p)}).val = 1;
    end
    tmp = fieldnames(hybrid.MF.params);
    idx = find(contains(tmp,'threshold'));
    for p = 1:length(idx)
        hybrid.MF.params.(tmp{idx(p)}).val = 1;
    end

    % Only select free parameters
    pnames = fieldnames(hybrid.params);
    np = length(pnames);
    freeidx = [];
    for p = 1:np
        if isnan(hybrid.params.(pnames{p}).val)
            freeidx = [freeidx p];
        end
    end
    pnames = pnames(freeidx);
    np = length(pnames);

    for mteach = {'mb-teach-mf','mf-teach-mb'}
        for mgain = {'new-old','old-new'}

            cc = cc+1;

            hybrid.teachtype = mteach{1};
            hybrid.gaintype = mgain{1};
        
            % save model to variable
            allmodels{cc} = hybrid;
        
            % Fit to subject data
            if ~onCluster
        
                this_nLL = nan(N,1);
                this_BIC = nan(N,1);
                this_ACC = nan(N,1);
        
                thisparams = nan(N,length(pnames));
                plotdata = nan(N,6,2);
                parfor s = 1:N
                
                    disp(['Subject ' subjects{s} ', MB = ' num2str(MB) ', MF = ' num2str(MF) '...'])
        
                    % Load data
                    d = parse_behav(subjects{s},dir_data);
                
                    % Fit parameters
                    [startvals,fitvals,err] = fit_startrange(hybrid,d,1,'both');
                
                    beststart = find(err==min(err),1,'first');
                    x = fitvals(:,beststart)';
                
                    % get predictions
                    [err, output, newd] = sim_model(d,hybrid,x,'both');
                    [aic,bic] = aicbic(err * (-1),length(x),sum(d.Forced==0));
                
                    % calculate accuracy
                    idx = d.Forced==0 & d.EV~=1;
                    acc = output.T.y(idx) == output.T.yhat(idx);
                
                    % log results
                    this_nLL(s,:) = err;
                    this_BIC(s,:) = bic;
                    this_ACC(s,:) = mean(acc);
                    thisparams(s,:) = x;
        
                    disp(['========= ' subjects{s} ' MB ' num2str(MB) ', MF ' num2str(MF) ' -->  accuracy = ' num2str(mean(acc))])
        
                    % plot
        %             %{
        %             figure
        %             subplot(1,2,1)
        %             y1 = nan(2,2); % choice, replaytype (rewarding, aversive)
        %             y1(1,1) = nanmean(newd.gain_rewarding(newd.Choice==1 & newd.Forced==0));
        %             y1(1,2) = nanmean(newd.gain_aversive(newd.Choice==1 & newd.Forced==0));
        %             y1(2,1) = nanmean(newd.gain_rewarding(newd.Choice==0 & newd.Forced==0));
        %             y1(2,2) = nanmean(newd.gain_aversive(newd.Choice==0 & newd.Forced==0));
        %             bar(y1')
        %             legend({'Rewarding','Aversive'})
        %             set(gca,'xticklabels',{'Approach','Avoid'})
        %             ylabel('Gain')
        %             subplot(1,2,2)
        %             y2 = nan(2,2); % choice, replaytype (rewarding, aversive)
        %             y2(1,1) = nanmean(newd.replayinitiated_rewarding(newd.Choice==1 & newd.Forced==0));
        %             y2(1,2) = nanmean(newd.replayinitiated_aversive(newd.Choice==1 & newd.Forced==0));
        %             y2(2,1) = nanmean(newd.replayinitiated_rewarding(newd.Choice==0 & newd.Forced==0));
        %             y2(2,2) = nanmean(newd.replayinitiated_aversive(newd.Choice==0 & newd.Forced==0));
        %             bar(y2')
        %             legend({'Rewarding','Aversive'})
        %             set(gca,'xticklabels',{'Approach','Avoid'})
        %             ylabel('Mean replay events')
        %             %}
                    
                    % Correlate with actual replay
                    thisreplay = replay(replay.Subject==str2double(subjects{s}),:);
                    cnames = {'gain_rewarding','gain_aversive','replayinitiated_rewarding','replayinitiated_aversive'};
                    y3 = nan(length(cnames),2);
                    for c = 1:length(cnames)
                        thisreplay.(cnames{c}) = nan(size(thisreplay,1),1);
                        for trl = 1:size(thisreplay,1)
                            idx = newd.ExpTrial==thisreplay.ExpTrial(trl);
                            thisreplay.(cnames{c})(trl) = newd.(cnames{c})(idx);
                        end
                        y3(c,1) = mean(thisreplay.(cnames{c})(thisreplay.Choice==1));
                        y3(c,2) = mean(thisreplay.(cnames{c})(thisreplay.Choice==0));
                    end
                    y3(5,1) = mean(thisreplay.Replay_rewarding(thisreplay.Choice==1));
                    y3(5,2) = mean(thisreplay.Replay_rewarding(thisreplay.Choice==0));
                    y3(6,1) = mean(thisreplay.Replay_aversive(thisreplay.Choice==1));
                    y3(6,2) = mean(thisreplay.Replay_aversive(thisreplay.Choice==0));
        
                    plotdata(s,:,:) = y3;
        
                    %{
                    figure
                    subplot(1,3,1)
                    bar(y3(5:6,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Actual replay')
                    subplot(1,3,2)
                    bar(y3(1:2,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Predicted gain')
                    subplot(1,3,3)
                    bar(y3(3:4,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Predicted replay events (average)')
                    %}
        
                end
        
                figure
                y = squeeze(nanmean(plotdata));
                subplot(1,3,1)
                bar(y(5:6,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Actual replay')
                subplot(1,3,2)
                bar(y(1:2,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Predicted gain')
                subplot(1,3,3)
                bar(y(3:4,:)'); legend({'Rewarding','Aversive'}); set(gca,'xticklabels',{'Approach','Avoid'}); title('Predicted replay events (average)')
                drawnow
        
                ACC(:,cc) = this_ACC';
                nLL(:,cc) = this_nLL';
                BIC(:,cc) = this_BIC';
        
                if length(allparams)<cc
                    allparams{cc} = struct();
                end
                for p = 1:length(pnames)
                    allparams{cc}.([pnames{p}])(s,:) = thisparams(:,p)';
                end
            else
        
                dir_save = fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',['cc' num2str(cc)]);
                if ~exist(dir_save)
                    mkdir(dir_save)
                end
        
                for s = 1:N
                    
                    filename = fullfile(dir_save,['modellingbatch-it' num2str(s) '.mat']);
                    model = hybrid;
                    info = [];
                    info.subject = subjects{s};
                    info.MB = MB;
                    info.MF = MF;
                    info.m = cc;
                    
                    thisd = parse_behav(subjects{s},dir_data);
        
                    save(filename,'model','thisd','numStart','info');
        
                end
        
                if onCluster
                    filename = fullfile(dir_save,['modellingbatch-it$SGE_TASK_ID.mat']);
                    generate_jobs_modelling(filename,'holly',cc);
                end
        
            end
        end
    end
end
save('D:\2020_RiskyReplay\results\modelling\hybridmodels.mat','allparams','allmodels','ACC','BIC','nLL');

% Group model performance
mnames = [];
midx = [];
for m = 1:length(allmodels)
    if strcmp(allmodels{m}.gaintype,'new-old') %&& strcmp(allmodels{m}.teachtype,'mf-teach-mb')
        midx = [midx m];
        mnames = [mnames {['MB' num2str(find(strcmp(extractfield(models,'name'),allmodels{m}.MB.name))),...
                        ', MF' num2str(find(strcmp(extractfield(models,'name'),allmodels{m}.MF.name))),...
                          ': '  allmodels{m}.teachtype(1:end-3) ', ' allmodels{m}.gaintype]}];
    end
end

[sorted,sortidx] = sort(sum(BIC(:,midx)));
figure
bar(sorted); set(gca,'xtick',1:length(midx)); set(gca,'xticklabels',mnames(sortidx)); ylabel('BIC'); title('Group sum')

figure
tmp = sum(BIC(:,midx));
bar([sum(tmp(contains(mnames,'mb-teach'))) sum(tmp(contains(mnames,'mf-teach')))]); set(gca,'xticklabels',{'mb-teach-mf','mf-teach-mb'}); ylabel('BIC'); title('')

% Individual model performance
[~,winner] = min(BIC(:,midx),[],2);
figure
histogram(winner)
set(gca,'xticklabels',mnames)
ylabel('Count'); title('Winning model per subject')

% Parameter estimates
cc = 18;
pnames = fieldnames(allparams{cc});
np = length(pnames);
pmatrix = nan(np,np);

figure
for p = 1:np

    y = allparams{cc}.(pnames{p});
    m = mean(y);
    sem = std(y)/sqrt(length(y));
    upper = m+sem;
    lower = m-sem;
    jitterx = repmat(p,length(y),1) + (rand(length(y),1)-0.5)*0.25;

    scatter(jitterx,y,'markerfacecolor','k','markeredgecolor','none'); hold on
    plot([p p],[upper lower],'r'); hold on
    scatter(p,m,80,'markerfacecolor','r'); hold on

    for p2 = 1:np
        pmatrix(p,p2) = corr(y',allparams{cc}.(pnames{p2})');
    end    
end
set(gca,'xtick',1:length(pnames))
set(gca,'xticklabels',pnames)

figure
imagesc(pmatrix);
set(gca,'xtick',1:length(pnames))
set(gca,'xticklabels',pnames)
set(gca,'ytick',1:length(pnames))
set(gca,'yticklabels',pnames)
caxis([-1 1])
colormap(colours(256,'redblue'))

% See whether model BIC or predicted replay events can predict actual replay
d = cell(1,N);
stackedreplay = [];
for m = 1:size(BIC,2)
    disp(['===== MODEL ' num2str(m) ' OF ' num2str(size(BIC,2)) '====='])
    for s = 1:N
        idx = replay.Subject==str2double(subjects{s});
        if any(idx)

            tmp = replay(idx,:);
            nRows = size(tmp,1);
            tmp.modelNum = repmat(m,nRows,1);
            tmp.mbNum = repmat(find(ismember(replace(allmodels{m}.MB.name,'_',' '),mnames)),nRows,1);
            tmp.mfNum = repmat(find(ismember(replace(allmodels{m}.MF.name,'_',' '),mnames)),nRows,1);
            tmp.BIC = repmat(BIC(s,m),nRows,1); % make it so that higher BIC is better fit
            tmp.teachtype = repmat({hybrid.teachtype},nRows,1);
            tmp.gaintype = repmat({hybrid.gaintype},nRows,1);
            
            pnames = fieldnames(allparams{m});
            params = [];
            for p = 1:length(pnames)
                params = [params allparams{m}.(pnames{p})(s)];
                if strcmp(pnames{p},'forgetting')
                    tmp.(['hybrid_MB_forgetting']) = repmat(allparams{m}.(pnames{p})(s),nRows,1); % duplicate it for now
                    tmp.(['hybrid_MF_forgetting']) = repmat(allparams{m}.(pnames{p})(s),nRows,1);
                else
                    tmp.(['hybrid_' pnames{p}]) = repmat(allparams{m}.(pnames{p})(s),nRows,1);
                end
            end
            
            if isempty(d{s})
                d{s} = parse_behav(subjects{s},dir_data);
            end
            [~, ~, newd] = sim_model(d{s},allmodels{m},params,'both');
            
            cnames = {'replayinitiated_rewarding','replayinitiated_aversive'};
            for c = cnames
                tmp.(c{1}) = nan(nRows,1);
                for trl = 1:nRows
                    tmp.(c{1})(trl) = newd.(c{1})(newd.ExpTrial==tmp.ExpTrial(trl));
                end
            end

            stackedreplay = [stackedreplay; tmp];
        end
    end
end
clear d;

df = [stackedreplay; stackedreplay; stackedreplay]; % rewarding, aversive, differential
df.Replay = [stackedreplay.Replay_rewarding; stackedreplay.Replay_aversive; stackedreplay.Replay_differential];
df.replayinitiated = [stackedreplay.replayinitiated_rewarding; stackedreplay.replayinitiated_aversive; stackedreplay.replayinitiated_rewarding-stackedreplay.replayinitiated_aversive];
df.PathType = [repmat({'rewarding'},size(stackedreplay,1),1); repmat({'aversive'},size(stackedreplay,1),1); repmat({'differential'},size(stackedreplay,1),1)];

df.modelNum = categorical(df.modelNum);

df = df(df.Lag>=20 & df.Lag<=90 & ~contains(df.PathType,'differential'),:);
% lme = fitglme(df, 'Replay~BIC*replayinitiated+(1|Subject:Lag)')

writetable(df,'hybridmodelswithreplay.csv');

% % Compare this to non-hybrid model performance
% 
% 
% % See whether either best group or best individual model's replay predictions correlate with actual replay
% replay = readtable('D:\2020_RiskyReplay\results\replay\replay_lme_planning.csv');
% 
% cc = 3;
% pnames = fieldnames(allparams{cc});
% np = length(pnames);
% 
% replay.replayinitiated_1 = nan(size(replay,1),1);
% replay.replayinitiated_2 = nan(size(replay,1),1);
% replay.gain_1 = nan(size(replay,1),1);
% replay.gain_2 = nan(size(replay,1),1);
% replay.Q_1 = nan(size(replay,1),1);
% replay.Q_2 = nan(size(replay,1),1);
% replay.QH_1 = nan(size(replay,1),1);
% replay.QH_2 = nan(size(replay,1),1);
% for s = 1:N
% 
%     % Load subject data
%     d = parse_behav(subjects{s},dir_data);
% 
%     % Get optimised parameters
%     params = nan(1,np);
%     for p = 1:np
%         params(p) = allparams{cc}.(pnames{p})(s);
%     end
% 
%     if s==1
%         for p = 1:np
%             replay.(['hybrid_' pnames{p}]) = nan(size(replay,1),1);
%         end
%     end
% 
%     idx = replay.Subject==str2double(subjects{s});
%     for p = 1:np
%         replay.(['hybrid_' pnames{p}])(idx) = repmat(params(p),sum(idx),1);
%     end
% 
%     % Fit model
%     [~, ~, newd] = sim_model(d,allmodels{cc},params);
% 
%     trials = unique(replay.ExpTrial(replay.Subject==str2double(subjects{s})));
%     for trl = 1:length(trials)
%         
%         idx1 = replay.Subject==str2double(subjects{s}) & replay.ExpTrial==trials(trl);
%         idx2 = contains(cellstr(newd.Subject),subjects{s}) & newd.ExpTrial==trials(trl);
% 
%         replay.replayinitiated_1(idx1) = newd.replayinitiated_1(idx2);
%         replay.replayinitiated_2(idx1) = newd.replayinitiated_2(idx2);
%         replay.gain_1(idx1) = newd.gain_1(idx2);
%         replay.gain_2(idx1) = newd.gain_2(idx2);
%         replay.Q_1(idx1) = newd.Q_1(idx2);
%         replay.Q_2(idx1) = newd.Q_2(idx2);
%         replay.QH_1(idx1) = newd.QH_1(idx2);
%         replay.QH_2(idx1) = newd.QH_2(idx2);
% 
%     end
% end
% 
% % GLM
% replay = replay(replay.Lag >= 20 & replay.Lag <= 90 & replay.RT>=5 & replay.RT<=30,:);
% 
% replay.replayinitiated_rewarding = nan(size(replay,1),1);
% idx = replay.nV_1 > replay.nV_2; replay.replayinitiated_rewarding(idx) = replay.replayinitiated_1(idx);
% idx = replay.nV_1 < replay.nV_2; replay.replayinitiated_rewarding(idx) = replay.replayinitiated_2(idx);
% 
% replay.replayinitiated_aversive = nan(size(replay,1),1);
% idx = replay.nV_1 < replay.nV_2; replay.replayinitiated_aversive(idx) = replay.replayinitiated_1(idx);
% idx = replay.nV_1 > replay.nV_2; replay.replayinitiated_aversive(idx) = replay.replayinitiated_2(idx);
% 
% replay.replayinitiated_cf = nan(size(replay,1),1);
% idx = replay.Choice==1; replay.replayinitiated_cf(idx) = replay.replayinitiated_aversive(idx);
% idx = replay.Choice==2; replay.replayinitiated_cf(idx) = replay.replayinitiated_rewarding(idx);
% 
% replay.Replay_counterfactual = nan(size(replay,1),1);
% idx = replay.Choice==1; replay.Replay_counterfactual(idx) = replay.Replay_differential(idx) * (-1);
% idx = replay.Choice==2; replay.Replay_counterfactual(idx) = replay.Replay_differential(idx);
% 
% replay.Replay_counterfactual = replay.Replay_counterfactual - nanmean(replay.Replay_counterfactual);
% replay.hybrid_replayalpha = replay.hybrid_replayalpha - nanmean(replay.hybrid_replayalpha);
% replay.hybrid_MB_forgetting = replay.hybrid_MB_forgetting - nanmean(replay.hybrid_MB_forgetting);
% replay.hybrid_MF_forgetting = replay.hybrid_MF_forgetting - nanmean(replay.hybrid_MF_forgetting);
% 
% lme = fitglme(replay, 'Replay_counterfactual~replayinitiated_cf+hybrid_replayalpha+hybrid_MF_forgetting+(1|Subject:Lag)')
% 
% writetable(replay,'D:\2020_RiskyReplay\results\modelling\hybridmodelswithreplay.csv');
% 
% % Plot replay predictions
% replay.gain_rewarding = nan(size(replay,1),1);
% idx = replay.nV_1 > replay.nV_2; replay.gain_rewarding(idx) = replay.gain_1(idx);
% idx = replay.nV_1 < replay.nV_2; replay.gain_rewarding(idx) = replay.gain_2(idx);
% 
% replay.gain_aversive = nan(size(replay,1),1);
% idx = replay.nV_1 < replay.nV_2; replay.gain_aversive(idx) = replay.gain_1(idx);
% idx = replay.nV_1 > replay.nV_2; replay.gain_aversive(idx) = replay.gain_2(idx);
% 
% replay.gain_counterfactual = nan(size(replay,1),1);
% replay.gain_counterfactual(replay.Choice==1) = replay.gain_aversive(replay.Choice==1);
% replay.gain_counterfactual(replay.Choice==2) = replay.gain_rewarding(replay.Choice==2);
% 
% figure
% for e = 1:2
% 
%     if e==1
%         evidx = replay.EV<0;
%     elseif e==2
%         evidx = replay.EV>0;
%     end
% 
%     y = nan(N,4);
%     for s = 1:N
%         idx = replay.Subject==str2double(subjects{s});
%         y(s,1) = mean(replay.gain_rewarding(idx & evidx & replay.Choice==1));
%         y(s,2) = mean(replay.gain_aversive(idx & evidx & replay.Choice==1));
%         y(s,3) = mean(replay.gain_rewarding(idx & evidx & replay.Choice==0));
%         y(s,4) = mean(replay.gain_aversive(idx & evidx & replay.Choice==0));
%     end
%     
%     subplot(1,2,e)
%     bar([nanmean(y(:,1:2)); nanmean(y(:,3:4))])
%     set(gca,'xticklabels',{'Approach','Avoid'})
%     legend({'Rewarding','Aversive'})
%     ylabel('Predicted replay')
% 
%     if e==1
%         title('Negative EV')
%     elseif e==2
%         title('Positive EV')
%     end
% 
% end
% 
% figure
% for cc = 1:size(modelcombos,1)
% 
%     y = nan(N,4);
%     for s = 1:N
%         idx = replay.Subject==str2double(subjects{s});
%         y(s,1) = mean(replay.gain_rewarding(idx & cidx & replay.Choice==1));
%         y(s,2) = mean(replay.gain_aversive(idx & cidx & replay.Choice==1));
%         y(s,3) = mean(replay.gain_rewarding(idx & cidx & replay.Choice==0));
%         y(s,4) = mean(replay.gain_aversive(idx & cidx & replay.Choice==0));
%     end
%     
%     subplot(1,3,c)
%     bar([nanmean(y(:,1:2)); nanmean(y(:,3:4))])
%     set(gca,'xticklabels',{'Approach','Avoid'})
%     legend({'Rewarding','Aversive'})
%     ylabel('Predicted replay')
% 
%     if c==1
%         title('Very certain')
%     elseif c==2
%         title('Moderately certain')
%     elseif c==3
%         title('Uncertain')
%     end
% 
% end
