function sequenceness = ss_build(Y,U,in)
% Generate sequenceness from timeseries of classifier predictions

%% Parameters

tp = [1 2; 2 3; 4 5; 5 6]; % transition pairs
nTransitions = size(tp,1);

nPerms = size(U,1);
nStates = 6;
nSamples = size(Y,1);

maxLag = 0.6; % in seconds
maxLag = maxLag / (1/in.Fs); % convert to samples

bins = 10; % consider 10 Hz alpha (or set to 'maxLag' to ignore)

%% Build

sequenceness = nan(size(U,1), 3, nTransitions, maxLag); % 1) fwd, bwd, diff, 3) transitions, 3) lag
for perm = 1:nPerms
    
    % Get this transition sequence
    thisPerm = U(perm,:);
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
    
%     % add path-to-path transition
%     M = zeros(nStates,nStates);
%     M(thisPerm(3),thisPerm(4)) = 1;
%     dM = [dM squash(M) squash(M')];

%     % account for differences in classifier accuracy
%     FC = zeros(nStates,nStates);
%     TC = zeros(nStates,nStates);
%     for t = 1:nTransitions % Transition
%         tp1 = thisPerm(tp(t,1));
%         tp2 = thisPerm(tp(t,2));
%         FC(tp1,tp2) = in.CA(tp1);
%         TC(tp1,tp2) = in.CA(tp2);
%         FC(tp2,tp1) = in.CA(tp2);
%         TC(tp2,tp1) = in.CA(tp1);
%     end
%     dM = [dM FC(:) TC(:)];
    
    % add autocorrelation and constant
    dM = [dM squash(eye(nStates)) squash(ones(nStates))];

    % Second-level GLM
    SEQ = pinv(dM) * betas';

    % Save to variable
    sequenceness(perm,1,:,:) = SEQ(1:nTransitions,:);
    sequenceness(perm,2,:,:) = SEQ(nTransitions+1:nTransitions*2,:);
    sequenceness(perm,3,:,:) = SEQ(1:nTransitions,:) - SEQ(nTransitions+1:nTransitions*2,:);
    
end


end