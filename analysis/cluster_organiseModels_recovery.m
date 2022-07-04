% function cluster_organiseModels_recovery()
% Gather all the output from the modelling

% addpath(genpath('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/'))

models = set_models();
nModel = length(models);
mnames = extractfield(models,'name');

nIterations = 50;

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

nLL = nan(nModels,nModels,nIterations);
AIC = nan(nModels,nModels,nIterations);
BIC = nan(nModels,nModels,nIterations);
recovery = cell(1,nModels);
for m1 = 1:nModels

    % load files
    filelist = dir(fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',...
        ['m' num2str(m1)],'modellingbatch*.mat'));

    pnames = fieldnames(models(m1).params);
    freeidx = [];
    for p = 1:length(pnames)
        if isnan(models(m1).params.(pnames{p}).val)
            freeidx = [freeidx p];
        end
    end
    np = length(freeidx);

    recovery{m1} = nan(nIterations,np,2);

    oldparams = cell(1,length(filelist));
    newparams = cell(1,length(filelist));
    thisaic = nan(length(filelist),1);
    thisbic = nan(length(filelist),1);
    thisnll = nan(length(filelist),1);
    allmnames = cell(length(filelist),1);
    paramidx = nan(length(filelist),3);
    parfor f = 1:length(filelist)

        disp(['M' num2str(m1) ': file ' num2str(f) ' of ' num2str(length(filelist))])

        tmp = load(fullfile(filelist(f).folder,filelist(f).name));
        try
            d = tmp.d;
        catch
            d = tmp.thisd;
        end
        idx = d.Forced==0;

        if m1 ~= tmp.info.m1
            error('Mismatching model number')
        end

        while ~isfield(tmp,'err')

            try
                warning(['--- model ' num2str(m1) ' -- ' filelist(f).name ' -- did not run; fitting it now...'])
                warning off; [~,tmp.fitvals,~] = fit_startrange(tmp.model,d,1,'one'); warning on;
    
                pnames = fieldnames(tmp.model.params);
                freeidx = [];
                for p = 1:length(pnames)
                    if isnan(tmp.model.params.(pnames{p}).val)
                        freeidx = [freeidx p];
                    end
                end
    
                tmp.fitvals = array2table(tmp.fitvals','variablenames',pnames(freeidx));
                [tmp.err, ~] = sim_model(d,tmp.model,table2array(tmp.fitvals),'one');
            catch
                disp(['((( some sort of error in while loop )))'])
            end

        end

        % Organise info
        it      = tmp.info.it;
        m2      = tmp.info.m2;
        err     = tmp.err;
        fitvals = table2array(tmp.fitvals);
        model   = tmp.model;

        % compute model fit
        beststart = find(err==min(err));
        x = mean(fitvals(beststart,:),1);

        [~,output,newd] = sim_model(d,model,x);

        if isempty(output)

            warning(['--- model ' num2str(m1) ' -- ' filelist(f).name ' -- did not run; fitting it now...'])
            [startvals,x,nLL] = fit_startrange(tmp.model,d,1,'both');
            [~,output,newd] = sim_model(d,model,x);
        end

        % Check that parameters of generative model don't make it indistinguishable from others
        pnames = fieldnames(models(m1).params)';
        freeidx = [];
        for p = 1:length(pnames)
            if isnan(models(m1).params.(pnames{p}).val)
                freeidx = [freeidx p];
            end
        end
        pnames = pnames(freeidx);

        oldparams{f} = tmp.info.origParams;

        % add to structures
        [aic,bic] = aicbic(output.nLL*(-1), length(x), size(output.T,1));
        nLL = output.nLL;

        thisaic(f) = aic;
        thisbic(f) = bic;
        thisnll(f) = nLL;

        newparams{f} = x;
        paramidx(f,:) = [m1 m2 it];
    end

    % Rearrange
    M2 = unique(paramidx(:,2))'; M2 = M2(~isnan(M2));
    iterations = unique(paramidx(:,3))'; iterations = iterations(~isnan(iterations));
    for m2 = M2
        for it = iterations
            idx = paramidx(:,1)==m1 & paramidx(:,2)==m2 & paramidx(:,3)==it;
            if any(idx)
                AIC(m1,m2,it) = thisaic(idx);
                BIC(m1,m2,it) = thisbic(idx);
                nLL(m1,m2,it) = thisnll(idx);
            end
        end
    end

    for it = iterations
        idx = paramidx(:,1)==m1 & paramidx(:,2)==m1 & paramidx(:,3)==it;
        if any(idx)
            recovery{m1}(it,:,1) = cell2mat(oldparams(idx)');
            recovery{m1}(it,:,2) = cell2mat(newparams(idx)');
        end
    end

    % Print results
    disp(['Model ' num2str(m1)])

    thisbic = squeeze(BIC(m1,:,:)); thisbic = thisbic(:,~any(isnan(thisbic)));
    [~,idx] = nanmin(thisbic);
    disp(['--- mean winning model = ' num2str(mean(idx==m1))])

    if mean(idx==m1) < 0.3
        warning('Low model recovery')
    end

    disp(['--- parameter recovery:'])
    y1 = squeeze(recovery{m1}(:,:,1)); y1 = y1(~all(isnan(y1),2),:);
    y2 = squeeze(recovery{m1}(:,:,2)); y2 = y2(~all(isnan(y2),2),:);
    corr(y1,y2)

end

save(fullfile('D:\2020_RiskyReplay\results\modelling','modelrecovery_evendistribution.mat'),...
    'recovery','AIC','BIC','nLL','mnames');

% end
