function [err, output, d] = sim_model(d,model,vec,randtype)
% [err, model] = sim_model(d,model,vec)
% 
% d = subject data table
% model = model specification (single model from output of set_models --> e.g., model = models(1))
% vec = vector of parameter values (in order specified in 'model' structure)
% randtype = 'one' = randomly pick paths (use for simulations), 'both' = take average of both (use for fitting)

if nargin <= 3
    randtype = 'both'; % 'one' or 'both'
end

if strcmp(model.pathtype,'random') && strcmp(randtype,'both')
    randPaths = 1:2;
else
    randPaths = 1;
end

%% Check input

% Get number of trials in data structure
nTrls = size(d,1);

% Make it so that 1 = approach, 0 = avoid
try
    if any(d.Choice==2)
        d.Choice(d.Choice==2) = 0;
    end
catch
end

% Check if parameters were provided (if not, use defaults)
if nargin==2

    params = [];

    params.alpha = 0.1;
    params.tau = 1;
    params.threshold = 1;
    params.rewthreshold = 1;
    params.lossthreshold = 1;
    params.replayalpha = 0.1;
    params.gain = 1;
    params.MB_forgetting = 0.1;
    params.MF_forgetting = 0.1;

else
    
    vec = vec(:)'; % make it a row

    % Cycle through parameters and note those with a free value
    pnames = fieldnames(model.params);
    np = length(pnames);

    freeidx = zeros(1,length(pnames));
    for p = 1:np
        freeidx(p) = isnan(model.params.(pnames{p}).val);
    end
    
    pnames = pnames(freeidx==1);
    np = length(pnames);

    % Create parameter fields and check values are within specified ranges
    for p = 1:np
        params.(pnames{p}) = vec(p);
        switch pnames{p}
            case {'alpha','replayalpha','MB_alpha','MF_alpha'}
                if vec(p) < 0 || vec(p) > 1
                    error('Learning rate (alpha) is out of bounds')
                end
            case {'tau','MF_tau','MB_tau'}
                if vec(p) <= 0
                    error('Inverse temperature (tau) is out of bounds')
                end
        end
    end

    % Get fixed parameters
    pnames = fieldnames(model.params);
    pnames = pnames(freeidx==0);
    np = length(pnames);

    % Create parameter fields and check values are within specified ranges
    for p = 1:np
        params.(pnames{p}) = model.params.(pnames{p}).val;
    end

end

% Check if choice data exists or not
existingdata = false;
colnames = d.Properties.VariableNames;
if ~any(contains(colnames,'Choice')) || ~any(contains(colnames,'Transition')) || ~any(contains(colnames,'Outcome'))
    
    d.Choice = nan(size(d,1),1);
    d.Transition = nan(size(d,1),1);
    d.Outcome = nan(size(d,1),1);

    % Put in equal number of transitions for forced choices at start of blocks
    blocks = unique([d.Practice d.Block],'rows');
    for b = 1:length(blocks)
        idx = d.Practice==blocks(b,1) & d.Block==blocks(b,2) & d.Forced~=0;
        d.Choice(idx) = 1;
        d.Transition(idx) = datasample([ones(1,sum(idx)/2) ones(1,sum(idx)/2)*2],sum(idx),'replace',false);
        d.Outcome(idx & d.Transition==1) = d.nV_1(idx & d.Transition==1);
        d.Outcome(idx & d.Transition==2) = d.nV_2(idx & d.Transition==2);
    end

else
    existingdata = true;
end

%% Estimate choices using model parameters

% Note neg positions per path, to match 'nCombo' in data table
nCombos = fliplr(combvec(1:3,1:3)');
nNegPos = size(nCombos,1);

% Note all possible probability conditions
pConditions = [.1 .3 .5 .7 .9];
nProb = length(pConditions);

% Combine probability conditions with neg combos
probNegCombos = combvec(1:nNegPos,1:nProb)';
nProbNeg = size(probNegCombos,1);

% Initialise variables
d.pV_1 = nan(size(d,1),1); % perceived path 1 value
d.pV_2 = nan(size(d,1),1); % perceived path 2 value
if contains(model.name,'random') && strcmp(model.pathtype,'random')
    d.randPathChoice = nan(size(d,1),1);
end

if strcmp(model.plantype,'hybrid')
    d.gain_1 = nan(size(d,1),1);
    d.gain_2 = nan(size(d,1),1);
    d.replayinitiated_1 = nan(size(d,1),1);
    d.replayinitiated_2 = nan(size(d,1),1);
    d.Q_1 = nan(size(d,1),1);
    d.Q_1 = nan(size(d,1),1);
    d.QH_1 = nan(size(d,1),1);
    d.QH_1 = nan(size(d,1),1);
end

V = [];
Q = []; 
switch model.calctype
    case 'calculation'
        if any(ismember(model.pathtype,{'rewarding','aversive'}))
            Q = [0 0];
        end
    case 'experience'
        switch model.qtype
            case 'perpath'
                Q = [0 0];
            case 'perneg'
                Q = zeros(3,2);
        end
    case 'actions'
        switch model.qtype
            case 'peraction'
                Q = [0 0];
            case 'perprob'
                Q = zeros(2,5); % action, probability type (10-90, 30-70, 50-50, 70-30, 90-10)
            case 'perneg'
                Q = zeros(2,nNegPos); % action, neg position
            case 'pernegprob'
                Q = zeros(2,5,nNegPos);
        end
end

pchoice = nan(nTrls,2,length(randPaths)); % cols: approach, avoid
errorcatch = false;
for RP = randPaths % do this twice - one for path 1 and one for path 2 - if randomly choosing paths (taking average probability)

    statevalues = zeros(2,3); 

    if length(randPaths)>1
        model.pathtype = ['path' num2str(randPaths(RP))];
    end

    for trl = 1:nTrls
    
        % Path appraisal
        prevstatevalues = statevalues;
        switch model.calctype
            case 'calculation'
                V = zeros(1,2);
                for path = 1:2
                    for st = 1:3
                        V(path) = round(V(path) + statevalues(path,st));
                        if mod(V(path),2)~=0 && st==nCombos(d.nCombo(trl),path)
                            V(path) = (-1) * V(path);
                        end
                    end
                end
        end
    
        % Path selection
        paths = [];
        if strcmp(model.plantype,'MB')
            switch model.pathtype
                case 'twopath'
                    paths = [1 2];
                case 'path1'
                    paths = 1;
                case 'path2'
                    paths = 2;
                case 'random'
                    switch randtype
                        case 'one'
                            paths = randi(2);
                            d.randPathChoice(trl) = paths;
                        case 'both'
                            paths = [1 2];
                    end
                case 'rewarding'
                    maxv = find(sum(Q,1)==max(sum(Q,1)));
                    if length(maxv)==1
                        paths = maxv;
                    else
                        paths = randi(2); % maxv(randi(2)); % if both paths are equivalent, either calculate BOTH (more accurate) or randomly pick ONE (less accurate)
                    end
                case 'aversive'
                    minv = find(sum(Q,1)==min(sum(Q,1)));
                    if length(minv)==1
                        paths = minv;
                    else
                        paths = randi(2); % minv(randi(2)); % if both paths are equivalent, either calculate BOTH (more accurate) or randomly pick ONE (less accurate)
                    end
            end
        end
    
        % Calculate expected value of approaching using perceived path value
        if strcmp(model.plantype,'MB')
    
            P = [d.P(trl) 1-d.P(trl)];
            switch model.calctype
                case 'calculation'
                    EVapp = V.*P;
                otherwise
                    if size(Q,1)==1
                        EVapp = Q(1:2).*P;
                    else
                        EVapp = [Q(nCombos(d.nCombo(trl),1),1) Q(nCombos(d.nCombo(trl),2),2)].*P;
                    end
            end

            d.pV_1(trl) = EVapp(1);
            d.pV_2(trl) = EVapp(2);
    
            EVapp = sum(EVapp(paths));

            % Add expected value of avoiding (i.e., threshold parameter)
            switch model.threshtype
                case 'single'
                    EV = [];
                    try
                        EV = [EVapp params.threshold];
                    catch
                        EV = [EVapp params.MB_threshold];
                    end
                case 'double'
                    EV = [];
                    try
                        if EVapp > 0
                            EV = [EVapp params.rewthreshold];
                        else
                            EV = [EVapp params.lossthreshold];
                        end
                    catch
                        if EVapp > 0
                            EV = [EVapp params.MB_rewthreshold];
                        else
                            EV = [EVapp params.MB_lossthreshold];
                        end
                    end
            end
        end
    
        % Determine choice
        switch model.plantype
            case 'MB'
                prob = exp(EV*params.tau) ./ sum(exp(EV*params.tau));
            case 'MF'
                thisQ = [];
                switch model.qtype
                    case 'peraction'
                        thisQ = Q;
                    case 'perprob'
                        thisQ = Q(:,pConditions==d.P(trl))';
                    case 'perneg'
                        thisQ = Q(:,d.nCombo(trl))';
                    case 'pernegprob'
                        thisQ = Q(:,pConditions==d.P(trl),d.nCombo(trl))';
                end
                prob = exp(thisQ*params.tau) ./ sum(exp(thisQ*params.tau));
            case 'bias'
                thisQ = [params.threshold 0]; % threshold is equivalent to a bias towards approaching
                prob = exp(thisQ*params.tau) ./ sum(exp(thisQ*params.tau));
        end
    
        if any(prob==1) || any(prob==0) || any(isnan(prob))
            errorcatch = true;
        end
    
        if d.Forced(trl)==0 % only make choices on free-choice trials
            % Insert into table
            pchoice(trl,:,RP) = prob;
        end
    
        % Update knowledge of path values
        if existingdata || d.Forced(trl)~=0
            transition = d.Transition(trl);
        else
            if d.Forced(trl)==0
                thischoice = datasample([1 0],1,'replace',false,'weights',pchoice(trl,:,RP)); % randomly choose action according to choice probabilities
                if thischoice==1
                    transition = datasample([1 2],1,'weights',[d.P(trl) 1-d.P(trl)]); % if approach, determine transition according to path probabilities
                    d.Outcome(trl) = d.(['nV_' num2str(transition)])(trl);
                else
                    transition = 0; % if avoid, no transition and guaranteed outcome of 1
                    d.Outcome(trl) = 1;
                end
                d.Transition(trl) = transition;
                d.Choice(trl) = thischoice;
            end
        end
   
        % If there was a transition, update perceived value based on outcomes
        if transition>0
    
            outcome = d.Outcome(trl);
    
            % Update values
            switch model.calctype
                case 'calculation'
                    if any(ismember(model.pathtype,{'rewarding','aversive'}))
                        Q(transition) = Q(transition) + params.alpha*(outcome - Q(transition));
                    end
                case 'experience'
                    switch model.qtype
                        case 'perpath'
                            Q(transition) = Q(transition) + params.alpha*(outcome - Q(transition));
                        case 'perneg'
                            Q(nCombos(d.nCombo(trl),transition),transition) = Q(nCombos(d.nCombo(trl),transition),transition) + params.alpha*(outcome - Q(nCombos(d.nCombo(trl),transition),transition));
                    end
            end
            if transition==1
                statevalues(transition,:) = [d.S1(trl) d.S2(trl) d.S3(trl)];
            elseif transition==2
                statevalues(transition,:) = [d.S4(trl) d.S5(trl) d.S6(trl)];
            end
        end
    
        % Update for either choice, if learning from actions
        switch model.calctype
            case 'actions'
    
                choice = [];
                if transition>0
                    choice = 1; % approach
                else
                    choice = 2; % avoid
                end
                outcome = d.Outcome(trl);
        
                switch model.qtype
                    case 'peraction'
                        Q(choice) = Q(choice) + params.alpha*(outcome - Q(choice));
                    case 'perprob'
                        Q(choice,pConditions==d.P(trl)) = Q(choice,pConditions==d.P(trl)) + params.alpha*(outcome - Q(choice,pConditions==d.P(trl)));
                    case 'perneg'
                        Q(choice,d.nCombo(trl)) = Q(choice,d.nCombo(trl)) + params.alpha*(outcome - Q(choice,d.nCombo(trl)));
                    case 'pernegprob'
                        Q(choice,pConditions==d.P(trl),d.nCombo(trl)) = Q(choice,pConditions==d.P(trl),d.nCombo(trl)) + params.alpha*(outcome - Q(choice,pConditions==d.P(trl),d.nCombo(trl)));
                end
        end
    
%         % Add forgetting
%         if ~isempty(Q)
%             Q = (1-params.forgetting)*Q;
%         end
% 
%         if ~strcmp(model.calctype,'calculation') && ~isempty(V)
%             V = (1-params.forgetting)*V;
%         else
%             prevstatevalues = (prevstatevalues + statevalues) / trl;
%             statevalues = (1-params.forgetting)*statevalues + (params.forgetting*prevstatevalues);
%         end
    end
end

%% Create output table (averaged over several iterations)

if errorcatch

    err = Inf;
    output = [];

else

    % To deal with stochastic models
    if length(randPaths)>1
        pchoice = squeeze(mean(pchoice,3)); % average over paths
    end

    % Ignore forced choice trials
    idx = d.Forced==0;
    
    T = array2table(nan(nTrls,5),'variablenames',{'y','yhat','prob','probApproach','probAvoid'});
    T.y = d.Choice;
    T.EV = d.EV;

    T.probApproach = pchoice(:,1);
    T.probAvoid = pchoice(:,2);
    
    T.prob = nan(nTrls,1);
    T.prob(T.y==1 & idx) = T.probApproach(T.y==1 & idx);
    T.prob(T.y==0 & idx) = T.probAvoid(T.y==0 & idx);
    
    T.yhat = nan(nTrls,1);
    for trl = 1:nTrls
        if idx(trl)
            try
                T.yhat(trl) = datasample([1 0],1,'replace',false,'weights',pchoice(trl,:));
            catch
                errorcatch = true;
            end
        end
    end

    err = -nansum(log(T.prob(d.Forced==0)));
    
    % Put into output structure
    output = [];
    output.T = T;
    output.nLL = err;

    if errorcatch
        err = Inf;
        output = [];
    end
    
%     disp(['Prediction accuracy = ' num2str(round(mean(d.Choice==T.yhat)*100,2)) '%'])
%     disp(['neg LL              = ' num2str(err)])
%     disp('-------------------------------------------')
end
end