% Gather all the output from the modelling

addpath('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/utils')

nModels = 28;
maxnp = 4;

nLL = [];
AIC = [];
BIC = [];
recovery = [];
paramlog = [];
for m1 = 1:nModels

    % load files
    filelist = dir(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/batch',...
        ['m' num2str(m1)],'modellingbatch*.mat'));

   for f = 1:length(filelist)

        disp(['M' num2str(m1) ': file ' num2str(f) ' of ' num2str(length(filelist))])

        tmp = load(fullfile(filelist(f).folder,filelist(f).name));

        if m1 ~= tmp.info.m1
            error('Mismatching model number')
        end

        it      = tmp.info.it;
        m2      = tmp.info.m2;
        err     = tmp.err;
        fitvals = tmp.fitvals;
        model   = tmp.model;
        d       = tmp.d;

        recovery{m1}(it,:,1) = tmp.info.origParams;
    
        % compute model fit
        beststart = find(err==min(err));
        x = mean(fitvals(:,beststart),2);

        if m1==m2
            recovery{m1}(it,:,2) = x;
        end
    
        [~,output] = sim_model(d,model,x);

        % add to structures
        [aic,bic] = aicbic(output.nLL*(-1), length(x), size(output.T,1));
        AIC(m1,m2,it) = aic;
        BIC(m1,m2,it) = bic;
        nLL(m1,m2,it) = output.nLL;
        paramlog(m1,m2,it,:) = [x' nan(1,maxnp-length(x))];

   end

    disp(['Model ' num2str(m1)])
    [~,idx] = min(squeeze(BIC(m1,:,:)));
    disp(['--- mean winning model = ' num2str(mean(idx==m1))])
    disp(['--- parameter recovery:'])
    corr(squeeze(recovery{m1}(:,:,1)),squeeze(recovery{m1}(:,:,2)))

end

save(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/behav/modelling','modelrecovery.mat'),...
    'recovery','AIC','BIC','nLL','paramlog');
