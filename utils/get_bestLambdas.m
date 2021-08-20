function [CV,lambdas] = get_bestLambdas(subject,trainTimes,thisnull)
%% Get best lambda for each classifier

dir_classifiers = 'D:\2020_RiskyReplay\data\meg\classifiers';

CV = [];
lambdas = [];
for t = 1:length(trainTimes)

    % Load cross-validation accuracy
    tmp = load(fullfile(dir_classifiers,subject,...
        ['cv_' subject '_t' num2str(trainTimes(t)) '_n' num2str(thisnull) '.mat']));
    thiscv = tmp.cv;

    % Average over folds
    thiscv = squeeze(nanmean(thiscv));

    % Get best lambda (averaged across conditions)
    [~,best_lambda] = max(nanmean(thiscv));
    lambdas(t) = best_lambda;

    % Get minimum accuracy across all 6 states
    CV(t,1) = min(thiscv(:,best_lambda)); 
    CV(t,2) = mean(thiscv(:,best_lambda)); 
    CV(t,3) = max(thiscv(:,best_lambda)); 

end

end