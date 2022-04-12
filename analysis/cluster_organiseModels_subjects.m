% Gather all the output from the modelling

addpath('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/utils')

% generate model structures
models = set_models();
mnames = extractfield(models,'name');
nModels = length(models);

numStart = 3;

subjects = {
    '012882'
    '013770'
    '015397'
    '018768'
    '027480'
    '088674'
    '097403'
    '147947'
    '220598'
    '263098'
    '383991'
    '391883'
    '396430'
    '503608'
    '506559'
    '521846'
    '663186'
    '680913'
    '706130'
    '707132'
    '795776'
    '832746'
    '909566'
    '945092'
    '957849'
    '989833'};
N = length(subjects);

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

for s = 1:N

    % load files
    filelist = dir(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/batch',...
        subjects{s},'modellingbatch*.mat'));

   for f = 1:length(filelist)

        disp([subjects{s} ': file ' num2str(f) ' of ' num2str(length(filelist))])

        tmp = load(fullfile(filelist(f).folder,filelist(f).name));

        m = find(strcmp(tmp.model.name,mnames));

        optim(m).startvals.start(s,:,:) = tmp.startvals;
        optim(m).startvals.nLL(s,:) = tmp.err;
        optim(m).startvals.fit(s,:,:) = tmp.fitvals;

        beststart = find(tmp.err==min(tmp.err));
        x = mean(tmp.fitvals(:,beststart),2);

        % predict data using optimised parameters
        [~,output] = sim_model(tmp.d,tmp.model,x);

        [aic,bic] = aicbic(-output.nLL,tmp.model.paraminfo.nFree,size(output.T,1));
            
        optim(m).nLL(s,:) = output.nLL;
        optim(m).params(s,:) = array2table(x');
        optim(m).aic(s) = aic;
        optim(m).bic(s) = bic;
        optim(m).acc(s) = mean(output.T.y(tmp.d.Forced==0)==output.T.yhat(tmp.d.Forced==0));

   end
end

save(fullfile('/data/holly-host/jmcfadyen/2020_RiskyReplay/behav/modelling','modelfitting.mat'),...
    'optim');
