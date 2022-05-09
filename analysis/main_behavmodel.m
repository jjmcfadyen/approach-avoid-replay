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

% generate model structures
models = set_models();
mnames = extractfield(models,'name');
nModels = length(models);

%% Simulate models to get performance

rng(1);

numStart = 1;
nIterations = 100;

optfit = [];
optfit.nLL = nan(nIterations,nModels);
optfit.acc = nan(nIterations,nModels);
optfit.params = nan(nIterations,nModels,4);
for it = 1:nIterations

    disp('====================================================================================')
    disp(['ITERATION ' num2str(it)])
    disp('====================================================================================')

    d = parse_behav(subjects{randi(N)},dir_data); % randomly select subject protocol for simulations
    
    d.Choice = nan(size(d,1),1); % make choices perfect
    d.Choice(d.EV>1 & d.Forced==0) = 1;
    d.Choice(d.EV<1 & d.Forced==0) = 0;
    d.Choice(d.EV==1 & d.Forced==0) = rand(sum(d.EV==1 & d.Forced==0))>.5;

    d.Transition(d.Choice==0) = 0; % update transitions
    d.Outcome(d.Choice==0) = 1;
    idx = find(d.Choice==1);
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
        
        % optimise parameters to perfect choice data
        [startvals,fitvals,err] = fit_startrange(models(m),d,numStart);

        beststart = find(err==min(err),1,'first');
        x = fitvals(:,beststart)';

        % get predictions
        [err, output] = sim_model(d,models(m),x);
        [aic,bic] = aicbic(err * (-1),models(m).paraminfo.nFree,sum(d.Forced==0));

        % calculate accuracy
        idx = d.Forced==0 & d.EV~=1;
        acc = output.T.y(idx) == output.T.yhat(idx);

        % log results
        nLL(m,:) = err;
        BIC(m,:) = bic;
        ACC(m,:) = mean(acc);

        if length(x) < size(optfit.params,3)
            x = [x nan(1,size(optfit.params,3)-length(x))];
        end
        params(m,:) = x;

    end

    optfit.nLL(it,:) = nLL;
    optfit.bic(it,:) = BIC;
    optfit.acc(it,:) = ACC;
    optfit.params(it,:,:) = params;

end

% Save
save(fullfile(dir_output,'comparestrategies.mat'),'optfit','models');

% Plot
figure
for i = 1:3

    subplot(1,3,i)

    if i==1
        y = optfit.nLL;
        thistitle = 'Negative Log Likelihood';
    elseif i==2
        y = optfit.bic;
        thistitle = 'BIC';
        y = y(:,end) - y;
    elseif i==3
        y = optfit.acc;
        thistitle = 'Accuracy';
    end

    m = mean(y);
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
    set(gca,'xticklabels',mnames)

end

%% Fit models

onCluster = true;
makePlots = false;
numStart = 3; % number of different starting values to use per parameter

% create empty structures for model optimisation output
if ~onCluster
    optim = struct();
    for m = 1:nModels
        np = models(m).paraminfo.nFree;
        optim(m).nLL = nan(N,1); % negative log likelihood
        optim(m).params = array2table(nan(N,np),'variablenames',models(m).paraminfo.names);
        optim(m).aic = nan(1,N);
        optim(m).bic = nan(1,N);
        optim(m).acc = nan(1,N);
        optim(m).startvals.start = nan(N,np,numStart^np); % subjects, parameters, starting combination index
        optim(m).startvals.fit = nan(N,np,numStart^np); % subjects, parameters, starting combination index
        optim(m).startvals.nLL = nan(N,numStart^np);
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

            [startvals,fitvals,err] = fit_startrange(models(m),d,numStart);
    
            optim(m).startvals.start(s,:,:) = startvals;
            optim(m).startvals.nLL(s,:) = err;
            optim(m).startvals.fit(s,:,:) = fitvals;
    
            beststart = find(err==min(err));
            x = mean(fitvals(:,beststart),2);
    
            % predict data using optimised parameters
            [~,output] = sim_model(d,models(m),x);
    
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
    
            [aic,bic] = aicbic(-output.nLL,models(m).paraminfo.nFree,size(output.T,1));
            
            optim(m).nLL(s,:) = output.nLL;
            optim(m).params(s,:) = array2table(x');
            optim(m).aic(s) = aic;
            optim(m).bic(s) = bic;
            optim(m).acc(s) = mean(output.T.y(idx)==output.T.yhat(idx));
        else
            
            cc = cc + 1;
            filename = fullfile(dir_save,['modellingbatch-it' num2str(cc) '.mat']);
            
            model = models(m);

            info = [];
            info.subject = subjects{s};
            info.m = m;

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
    save(fullfile(dir_output,'modelfitting.mat'),'models','optim');
end

%% Plot

modelorder = 1:nModels;

excludedSubjects = {'263098','680913'};
includedidx = ~ismember(subjects,excludedSubjects);

% get negative log likelihoods
nLL = nan(N,nModels);
for m = 1:nModels
    nLL(:,m) = optim(m).nLL;
end

% calculate AIC and BIC
AIC = nan(N,nModels);
BIC = nan(N,nModels);
for m = 1:nModels
    AIC(:,m) = optim(m).aic;
    BIC(:,m) = optim(m).bic;
end

% group model evidence
figure
bar(sum(BIC(includedidx,end)) - sum(BIC(includedidx,modelorder)))
title('group model evidence')
set(gca,'xtick',1:nModels)
set(gca,'xticklabels',mnames(modelorder))

% subject winning model histogram
[~,winner] = min(BIC(includedidx,modelorder),[],2);
figure
histogram(winner(includedidx))
title('Winning model per subject')
set(gca,'xtick',1:nModels)
set(gca,'xticklabels',mnames(modelorder))

% make table to put into replay analyses
winnertable = array2table(subjects,'variablenames',{'Subject'});
winnertable.bModelNum = winner;

winnertable.bCalculationType = cell(N,1);
winnertable.bPathChoice = cell(N,1);
for s = 1:N
    thismodel = mnames{winnertable.bModelNum(s)};
    if contains(thismodel,'Q') || strcmp(thismodel,'null')
        winnertable.bCalculationType{s} = 'caching';
    else
        winnertable.bCalculationType{s} = 'calculating';
    end
    if contains(thismodel,'optimal')
        winnertable.bPathChoice{s} = 'both';
    else
        winnertable.bPathChoice{s} = 'one';
    end
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

% show individual subjects, ordered by accuracy
[sorted,sortidx] = sort(acc,'descend');

figure
for s = 1:N
    subplot(4,7,s)
    bic = BIC(sortidx(s),:);
    bic = bic(strcmp(mnames,'null')) - bic(:);
    bic = bic(modelorder);
%     tmp = [[2:10]' bic(2:end)];
    bar(bic);
    title([num2str(round(acc(sortidx(s))*100)) '% - ' mnames{find(bic==max(bic))}])
    set(gca,'ticklength',[0 0])
end

% parameter estimates
figure
for m = 1:nModels

    subplot(4,8,m)
    
    paramnames = models(m).paraminfo.names;
    np = models(m).paraminfo.nFree;
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

onCluster = true;

nIterations = 50;
numStart = 1;

maxnp = 4;

% run procedure
if ~onCluster
    nLL = nan(nModels,nModels,nIterations);
    AIC = nan(nModels,nModels,nIterations);
    BIC = nan(nModels,nModels,nIterations);
    recovery = cell(1,nModels);
    paramlog = nan(nModels,nModels,nIterations,maxnp);
end

D = cell(1,N);
for s = 1:N
    d = parse_behav(subjects{randi(N)},dir_data); % randomly select subject protocol for simulations
    d = d(:,1:find(contains(d.Properties.VariableNames,'Choice'))-1);
    D{s} = d;
end

for m1 = 1:nModels

    disp('=====================================================')
    disp(['MODEL ' num2str(m1)])
    disp('=====================================================')
    pnames = models(m1).paraminfo.names;

    if onCluster
        batchdir = fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',['m' num2str(m1)]);
        if ~exist(batchdir)
            mkdir(batchdir);
        end
        cc = 0;
    end

    % Generate choices using lots of different sets of parameters
    origP = nan(nIterations,models(m1).paraminfo.nFree); % 1 = iteration number, 2 = parameter number
    for it = 1:nIterations
     
        % generate parameter set
        itparams = nan(models(m1).paraminfo.nFree,1);
        for p = 1:models(m1).paraminfo.nFree
            
            fitparams = table2array(optim(m1).params(:,p));
            thisparam = normrnd(mean(fitparams),std(fitparams));
            thisrange = [models(m1).params.(models(m1).paraminfo.names{p}).lb models(m1).params.(models(m1).paraminfo.names{p}).ub];
            while thisparam < thisrange(1) || thisparam > thisrange(2)
                thisparam = normrnd(mean(fitparams),std(fitparams)); % make sure the randomly generated number is within the correct range
            end

%             if isinf(models(m1).params.(pnames{p}).ub)
%                 thisparam = unifrnd(realmin,5);
%             else
%                 thisparam = unifrnd(models(m1).params.(pnames{p}).lb,models(m1).params.(pnames{p}).ub);
%             end

            itparams(p,1) = thisparam;
        end
        origP(it,:) = itparams;
    end

    if ~onCluster
        newP = nan(nIterations,models(m1).paraminfo.nFree); % 1 = iteration number, 2 = parameter number
        parfor it = 1:nIterations
    
            % simulate new data
            d = D{randi(N)} % randomly select subject protocol for simulations

            % generate predicted data
            [~, output, thisd] = sim_model(d,models(m1),squeeze(origP(it,:)));
    
            % predict data using each model
            for m2 = 1:nModels

                [startvals,fitvals,err] = fit_startrange(models(m2),thisd,numStart);
        
                beststart = find(err==min(err));
                x = mean(fitvals(:,beststart),2);
    
                if m2==m1
                    newP(it,:) = x;
                end
    
                [~,output] = sim_model(thisd,models(m2),x);
    
                [aic,bic] = aicbic(output.nLL*(-1), models(m2).paraminfo.nFree, size(output.T,1));
                AIC(m1,m2,it) = aic;
                BIC(m1,m2,it) = bic;
                nLL(m1,m2,it) = output.nLL;
                paramlog(m1,m2,it,:) = [x' nan(1,maxnp-length(x))];
            end
        end
    else

        for it = 1:nIterations

            % simulate new data
            d = D{randi(N)};% randomly select subject protocol for simulations

            % generate predicted data
            [~, output, thisd] = sim_model(d,models(m1),squeeze(origP(it,:)));

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

%% Compute metrics

modelorder = 1:nModels;

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
        y(m1,m2) = sum(idx==m2);
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
imagesc(y(modelorder,modelorder)')
caxis([0 1])
colormap('gray')
set(gca,'ticklength',[0 0])
set(gca,'xtick',1:length(modelorder))
set(gca,'xticklabels',mnames(modelorder))
set(gca,'ytick',1:length(modelorder))
set(gca,'yticklabels',mnames(modelorder))
xlabel('Generative model')
ylabel('Winning predictive model')

%% Get predicted path replay (based on counterfactual utility)

% Select subjects
excludedSubjects = {'263098','680913'};
includedidx = ~ismember(subjects,excludedSubjects);

thesesubjects = subjects(includedidx);
thisN = length(thesesubjects);

% Get winning models
midx = [1 10 19]; % all two-path models
BIC = nan(thisN,length(midx));
for m = 1:length(midx)
    BIC(:,m) = optim(midx(m)).bic(includedidx);
end

% Make table
counterfactual = [];
for s = 1:thisN

    % Get behavioural data
    d = parse_behav(thesesubjects{s},dir_data); 

    % Get best learning model
    thesemodels = midx([1 find(BIC(s,2:end) == min(BIC(s,2:end))) + 1]);

    % Get path utility per model
    d.utility_counter_mCalculate = nan(size(d,1),1);
    d.utility_counter_mLearn = nan(size(d,1),1);
    for m = thesemodels
        if m==1
            msubname = 'mCalculate';
        else
            msubname = 'mLearn';
        end
        [~, ~, thisd] = sim_model(d,models(m),table2array(optim(m).params(find(contains(subjects,thesesubjects{s})),:)));
        pathvals = [thisd.pV_1 thisd.pV_2];
        for trl = 1:size(d,1)
            
            if ~all(isnan(pathEV))
                if d.Choice(trl)==1
                    d.(['utility_counter_' msubname])(trl) = pathvals(pathvals==min(pathEV));
                else
                    d.(['utility_counter_' msubname])(trl) = pathvals(pathvals==max(pathEV));
                end
            end
        end
    end
    counterfactual = [counterfactual; d];
end

% Exclude forced-choice & slow/fast RTs
counterfactual = counterfactual(counterfactual.RT>= 5 & counterfactual.RT <= 30 & counterfactual.Forced==0,:);

% Save to table for analysis in R
writetable(counterfactual,'D:\2020_RiskyReplay\results\modelling\counterfactualtable.csv');
