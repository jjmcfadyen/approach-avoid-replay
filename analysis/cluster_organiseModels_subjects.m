% function cluster_organiseModels_subjects()
% % Gather all the output from the modelling
% 
% addpath('/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/utils')

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
for s = 1:N

    % load files
    filelist = dir(fullfile('D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch',...
        subjects{s},'modellingbatch*.mat'));

   for m = 1:length(filelist)

        disp([subjects{s} ': file ' num2str(m) ' of ' num2str(length(filelist))])

        tmp = load(fullfile(filelist(m).folder,filelist(m).name));

        if ~isfield(tmp,'startvals')
            warning(['Model did not run for s=' num2str(s) ', file = ' fullfile(filelist(m).folder,filelist(m).name)])
            fit_startrange(fullfile(filelist(m).folder,filelist(m).name));
            tmp = load(fullfile(filelist(m).folder,filelist(m).name));
        end

        if s==1
            optim(m).model = tmp.model;
        else
            if ~strcmp(optim(m).model.name,tmp.model.name)
                error('Different model names!')
            end
        end

        optim(m).startvals.start(s,:,:) = tmp.startvals;
        optim(m).startvals.nLL(s,:) = tmp.err;
        optim(m).startvals.fit(s,:,:) = table2array(tmp.fitvals);

        beststart = find(tmp.err==min(tmp.err));
        x = mean(table2array(tmp.fitvals(beststart,:)),1);

        % predict data using optimised parameters
%         if contains(tmp.model.name,'random')
%             disp(['Running random iterations for s' num2str(s) ', m' num2str(m) '...'])
%             nIterations = 100;
%             acc = nan(1,nIterations);
%             bic = nan(1,nIterations);
%             nLL = nan(1,nIterations);
%             params = nan(length(x),nIterations);
%             parfor it = 1:nIterations
%                 warning off
%                 [~,x,~] = fit_startrange(tmp.model,tmp.d,1,'one');
%                 warning on
%                 [~,output] = sim_model(tmp.d,tmp.model,x,'one');
%                 [A,B] = aicbic(-output.nLL,length(x),size(output.T,1));
%                 acc(it) = mean(output.T.y(tmp.d.Forced==0)==output.T.yhat(tmp.d.Forced==0));
%                 bic(it) = B;
%                 nLL(it) = output.nLL;
%                 params(:,it) = x;
%             end
%             acc = mean(acc);
%             bic = mean(bic);
%             nLL = sum(nLL);
%             x = mean(params,2)';
%         else
            [~,output] = sim_model(tmp.d,tmp.model,x,'both');
            [aic,bic] = aicbic(-output.nLL,length(x),size(output.T,1));
            acc = mean(output.T.y(tmp.d.Forced==0)==output.T.yhat(tmp.d.Forced==0));
            nLL = output.nLL;
%         end
            
        optim(m).nLL(s,:) = nLL;
        if s==1
            optim(m).params = array2table(nan(N,length(x)),'variablenames',tmp.fitvals.Properties.VariableNames);
        end
        optim(m).params(s,:) = array2table(x,'variablenames',optim(m).params.Properties.VariableNames);
        optim(m).aic(s) = aic;
        optim(m).bic(s) = bic;
        optim(m).acc(s) = acc;

   end
end

save(fullfile('D:\2020_RiskyReplay\results\modelling',['modelfitting_randtype-' tmp.info.randtype '.mat']),'optim');


% end