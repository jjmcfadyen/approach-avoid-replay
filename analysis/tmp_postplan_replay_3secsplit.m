planningType = 'post'; % 'during' planning (main analysis) or 'post' planning (transition to outcome)
temporalModulation = false;
dir_replay = 'D:\2020_RiskyReplay\data\meg\replay\withoutintercept_postplanning';

for s = 1:N
    
    dir_save = fullfile(dir_replay,subjects{s});
    if ~exist(dir_save)
        mkdir(dir_save);
    end
    
    % Get task data
    load(fullfile(dir_meg,['7_merged_ds-' num2str(Fs) 'Hz'],[subjects{s} '_task_' num2str(Fs) 'Hz.mat'])); % loads 'merged' variable
    
    if ~isfield(merged,'fsample')
        merged.fsample = Fs;
    end
    
    for e = 1:2 % first 3 seconds, last 3 seconds

        data = merged;
        for trl = 1:length(data.trial)
            if e==1
                idx = data.time{trl} >= 0 & data.time{trl} <= 3;
            elseif e==2
                idx = data.time{trl} >= data.time{trl}(end)-3;
            end
            data.trial{trl} = data.trial{trl}(:,idx);
            data.time{trl} = data.time{trl}(idx);
        end
        
        % Get replay for each classifier training time
        thesetimes = 110:10:140;
        for t = 1:length(thesetimes)

    %         % Create job for cluster
    %         generate_jobs_replay(subjects{s},thesetimes(t),thisnull,lambdas(s,trainTimes==thesetimes(t)));
    %         
    %         disp('=========================================================')
    %         disp(['=== ' subjects{s} ', TRAINING TIME ' num2str(t) ' of ' num2str(length(thesetimes)) ' ======================'])
    %         disp('=========================================================')
    %         
    %         tic
    %         
    % %         % build & save classifier
    % %         filename = fullfile(dir_classifiers,subjects{s},['data_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']);
    % %         classifier = build_classifier(filename);
    % %         save(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'classifier');
            load(fullfile(dir_classifiers,subjects{s},['classifier_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']));

            % get replay using best lambdas from classifier & save
            thislambda = lambdas(s,trainTimes==thesetimes(t));
            classifier.betas = squeeze(classifier.betas(:,:,thislambda));
            classifier.intercepts = classifier.intercepts(:,thislambda);
            if ~temporalModulation
                replay = compute_replay(data,classifier,U);
            else
                replay = compute_replay_tm(data,classifier,U); 
            end

            timetag = '';
            if temporalModulation
                timetag = '_time-modulated';
            end
            
            if e==1
                etag = '_transition';
            elseif e==2
                etag = '_outcome';
            end
            
            save(fullfile(dir_save,['replay' etag timetag '_' subjects{s} '_t' num2str(thesetimes(t)) '_n' num2str(thisnull) '.mat']),'replay');

    %         toc

        end
    end
end