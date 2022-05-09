function sequenceness = compute_replay(data,classifier,nullperms)

%% Get info

nPerms = size(nullperms,1);
nStates = 6;
nTrls = length(data.trial);

%% Settings

% Transitions
tp = [1 2; 2 3; 4 5; 5 6]; % transition pairs
nTransitions = size(tp,1);

% Lags
maxLag = 0.6; % in seconds
maxLag = maxLag / (1/data.fsample); % convert to samples

bins = data.fsample / 10; % consider 10 Hz alpha (or set to 'maxLag' to ignore)

%% Compute replay per trial

sequenceness = nan(nTrls,nPerms,3,nTransitions,maxLag);
for trl = 1:nTrls
   
    disp(['Computing replay for trial ' num2str(trl) ' of ' num2str(nTrls) '...'])
    
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
    Y = X*C.betas';
    
    thisreplay = [];
    for perm = 1:nPerms
        
        % Get this transition sequence
        thisPerm = nullperms(perm,:);
        if perm == 1 && sum(thisPerm == [1:nStates]) ~= nStates
            error('First permutation is not the main permutation')
        end

        % Create toeplitz matrix
        TM = nan(nSamples,nStates*maxLag);
        cc = linspace(0,size(TM,2),nStates+1);
        warning off
        for st = 1:nStates
            tmp = toeplitz(Y(:,st),zeros(maxLag+1,1));
            TM(:,cc(st)+1:cc(st+1)) = tmp(:,2:end);
        end
        warning on

        % First-level GLM
        betas = nan(nStates*maxLag, nStates);
        for ilag = 1:bins
            idx = (1:bins:nStates*maxLag) + ilag - 1;
            tmp = pinv([TM(:,idx) ones(size(TM,1),1)])*Y;
            betas(idx,:) = tmp(1:end-1,:);
        end

        betas = reshape(betas,[maxLag nStates^2]);

        % Design matrix - transitions in this permutation
        dM = [];
        for d = 1:2 % Direction: forward, backward
            for t = 1:nTransitions % Transition
                M = zeros(nStates,nStates);
                tp1 = thisPerm(tp(t,1));
                tp2 = thisPerm(tp(t,2));
                M(tp1,tp2) = 1;
                if d == 2
                    M = M'; % transpose for backward transition
                end
                dM = [dM, M(:)];
            end
        end
        
        % add autocorrelation and constant
        dM = [dM squash(eye(nStates)) squash(ones(nStates))];

        % Second-level GLM
        SEQ = pinv(dM) * betas';

        % Save to variable
        thisreplay(perm,1,:,:) = SEQ(1:nTransitions,:);
        thisreplay(perm,2,:,:) = SEQ(nTransitions+1:nTransitions*2,:);
        thisreplay(perm,3,:,:) = SEQ(1:nTransitions,:) - SEQ(nTransitions+1:nTransitions*2,:);
        
    end
    
    sequenceness(trl,:,:,:,:) = thisreplay;
    
end

disp('Finished!')

end