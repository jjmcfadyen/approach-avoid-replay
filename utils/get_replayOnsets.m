function [onsets seqevidence] = get_replayOnsets(data,classifier,lag)
% Identify replay onsets

%% Get info

nStates = 6;
nTrls = length(data.trial);

transitions = [1 2; 2 3; 4 5; 5 6];

%% Get onsets

onsets = [];
seqevidence = cell(1,nTrls);
for trl = 1:nTrls
   
    C = classifier;
    
    % Scale data
    X = data.trial{trl}'; % channels x samples
    X = X ./ prctile(abs(X(:)),95);
    
    if any(isnan(X(:)))
        idx = ~any(isnan(X));
        X = X(:,idx);
        C.betas = C.betas(:,idx);
    end

    nSamples = size(X,1);
    
    % Apply classifier to get predicted data
%     Y = normr(1 ./ (1 + exp(-(X*C.betas' + repmat(C.intercepts', [size(X,1) 1])))));
%     Y = X*C.betas';
    Y = normr(1 ./ (1 + exp(-(X*C.betas'))));
%     Y = normalise(X*C.betas');

    % Get time-shifted reactivation matrix
    threshold = quantile(Y(:),.95); % 95th percentile of overall reactivation matrix
    thisseqevidence = nan(nSamples,size(transitions,1));
    theseonsets = nan(nSamples,size(transitions,1));
    for t = 1:size(transitions,1)
        TM = nan(nSamples,length(lag)+1);
        TM(:,1) = Y(:,transitions(t,1)); % this state's reactivation
        for ilag = 1:length(lag)
            lagsamples = round(lag(ilag)/1000*data.fsample); % in samples
            TM(:,ilag+1) = [Y((lagsamples+1):end,transitions(t,2)); nan(lagsamples,1)]; % reactivation of next state at this lag
        end
        [bestamp,bestlag] = max(TM(:,2:end),[],2);
        theseonsets(:,t) = TM(:,1) > threshold & bestamp > threshold;
        thisseqevidence(:,t) = TM(:,1) + bestamp;
    end
    seqevidence{trl} = thisseqevidence;
    
    binonsets_overall = nansum(theseonsets,2) > 0;
    binonsets_perpath = [nansum(theseonsets(:,1:2),2) > 0, nansum(theseonsets(:,3:4),2) > 0];
    
    % remove any onsets where there are onsets in the previous 100 ms
    window = 100; % in ms
    for path = 1:3 % path 1, path 2, overall
        if path<3
            idx = find(binonsets_perpath(:,path));
        else
            idx = find(binonsets_overall);
        end
        for i = 1:length(idx)
            
            % subset data in preceding 100ms
            samplewindow = (idx(i)-(window/1000*data.fsample)):(idx(i)-1);
            samplewindow = samplewindow(samplewindow>0 & samplewindow<size(binonsets_perpath,1)); % remove negative samples
            
            if path<3
                if any(binonsets_perpath(samplewindow,path)==1)
                    binonsets_perpath(samplewindow,path) = 0;
                end
            else
                if any(binonsets_overall(samplewindow,:)==1)
                    binonsets_overall(samplewindow,:) = 0;
                end
            end
        end
    end
    
    % Add to table  
    nEvents = sum(binonsets_perpath(:)) + sum(binonsets_overall(:));
    
    thistable = [];
    if data.trialinfo(trl,1)==0
        thistable.Practice = repmat(1,nEvents,1);
        thistable.Block = repmat(1,nEvents,1);
    else
        thistable.Practice = repmat(0,nEvents,1);
        thistable.Block = repmat(data.trialinfo(trl,1),nEvents,1);
    end
    thistable.Trial = repmat(data.trialinfo(trl,2),nEvents,1);
    
    thistable.Path = [ repmat(0,sum(binonsets_overall),1); repmat(1,sum(binonsets_perpath(:,1)),1); repmat(2,sum(binonsets_perpath(:,2)),1)];
    thistable.Onset_sample = [find(binonsets_overall); find(binonsets_perpath(:,1)); find(binonsets_perpath(:,2))];
    thistable.Onset_time = data.time{trl}(thistable.Onset_sample)';
    
    thistable.Duration_sample = repmat(nSamples,nEvents,1); % total time in the trial window
    
    onsets = [onsets; sortrows(struct2table(thistable),'Onset_sample')];
    
end

end