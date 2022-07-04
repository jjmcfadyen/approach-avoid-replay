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

if ~strcmp(model.plantype,'hybrid')
    if strcmp(model.pathtype,'random') && strcmp(randtype,'both')
        randPaths = 1:2;
    else
        randPaths = 1;
    end
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

% Account for hybrid models
if strcmp(model.plantype,'hybrid')
    MF = model.MF;
    MB = model.MB;
else
    MF = model;
    MB = model;
end

% Initialise variables
d.pV_1 = nan(size(d,1),1); % perceived path 1 value
d.pV_2 = nan(size(d,1),1); % perceived path 2 value

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
if ~strcmp(model.plantype,'hybrid')
    switch model.calctype
        case 'calculation'
            if any(ismember(model.pathtype,{'rewarding','aversive'}))
                V = [0 0];
            end
        case 'experience'
            switch model.qtype
                case 'perpath'
                    V = [0 0];
                case 'perneg'
                    V = zeros(3,2);
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
else

    % Use 'Q' for model free values
    switch model.MF.calctype
        case 'calculation'
            error('It does not make sense to have a model free system that calculates values')
        case 'experience'
            switch model.MF.qtype
                case 'perpath'
                    Q = [0 0 0]; % path 1, path 2, safe outcome
                case 'perneg'
                    Q = zeros(3,3); % path (1, 2, safe), neg position (1 to 3)
                    Q(3,2:3) = NaN;
            end
        case 'actions'
            switch model.MF.qtype
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

    % Duplicate to get hybrid Q values
    QH = Q;

    % Use 'V' for model-based values
    switch model.MB.calctype
        case 'calculation'
            if any(ismember(model.MB.pathtype,{'rewarding','aversive'}))
                V = [0 0];
            end
        case 'actions'
            error('It does not make sense for the model-based system to learn action-state values')
    end
end

if strcmp(model.plantype,'hybrid')
    if strcmp(model.MB.name, model.MF.name)
        V = []; % use Q values for both
    end
end

pchoice = nan(nTrls,2,length(randPaths)); % cols: approach, avoid
errorcatch = false;
for RP = randPaths % do this twice - one for path 1 and one for path 2 - if randomly choosing paths (taking average probability)

    statevalues = zeros(2,3); 

    if length(randPaths)>1
        MB.pathtype = ['path' num2str(randPaths(RP))];
    end

    for trl = 1:nTrls
    
        % Path appraisal
        prevstatevalues = statevalues;
        switch MB.calctype
            case 'calculation'
                if ~any(ismember(MB.pathtype,{'rewarding','aversive'}))
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
        end
    
        % Path selection
        paths = [];
        if strcmp(MB.plantype,'MB')
            switch MB.pathtype
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
                        case 'both'
                            paths = [1 2];
                    end
                case 'rewarding'
                    maxv = find(sum(V,1)==max(sum(V,1)));
                    if length(maxv)==1
                        paths = maxv;
                    else
                        paths = [1 2]; % maxv(randi(2)); % if both paths are equivalent, either calculate BOTH (more accurate) or randomly pick ONE (less accurate)
                    end
                case 'aversive'
                    minv = find(sum(V,1)==min(sum(V,1)));
                    if length(minv)==1
                        paths = minv;
                    else
                        paths = [1 2]; % minv(randi(2)); % if both paths are equivalent, either calculate BOTH (more accurate) or randomly pick ONE (less accurate)
                    end
            end
        end
    
        % Calculate expected value of approaching using perceived path value
        if strcmp(MB.plantype,'MB')
    
            P = [d.P(trl) 1-d.P(trl)];
            if isempty(V)
                if size(Q,1)==1
                    EVapp = Q(1:2).*P;
                else
                    EVapp = [Q(1,nCombos(d.nCombo(trl),1)) Q(2,nCombos(d.nCombo(trl),2))].*P;
                end
            else
                EVapp = V.*P;
            end
    
            if ~strcmp(model.plantype,'hybrid')
                if strcmp(MB.calctype,'experience') && strcmp(MB.qtype,'perneg')
                    EVapp = [EVapp(nCombos(d.nCombo(trl),1),1) EVapp(nCombos(d.nCombo(trl),2),2)];
                end
        
    %             if strcmp(MB.pathtype,'random') && strcmp(randtype,'both') % need to consider both possibilities
    %                 EVapp = EVapp;
    %             else
                    EVapp = sum(EVapp(paths));
    %             end
    
                % Add expected value of avoiding (i.e., threshold parameter)
                switch MB.threshtype
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
            else
    
                EV = [];
                switch MB.threshtype
                    case 'single'
                        try
                            EV = [EVapp params.threshold];
                        catch
                            EV = [EVapp params.MB_threshold];
                        end
                    case 'double'
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
        end
    
        % Determine choice
        if ~strcmp(model.plantype,'hybrid')
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
        else
    
            % Get Q value for model-free system and compute "hybrid" Q value
            thisQ = [];
            switch MF.qtype
                case 'perpath' % technically a "model-based" model
                    thisQ = Q;
                case 'peraction'
                    thisQ = [Q(1)/2 Q(1)/2 Q(2)];
                case 'perprob'
                    thisQ = Q(:,pConditions==d.P(trl))';
                    thisQ = [thisQ(1) thisQ(1) thisQ(2)];
                case 'perneg'
                    thisQ = [Q(1,nCombos(d.nCombo(trl),1)) Q(2,nCombos(d.nCombo(trl),2)) Q(3,1)];
                case 'pernegprob'
                    thisQ = Q(:,pConditions==d.P(trl),d.nCombo(trl))';
                    thisQ = [thisQ(1) thisQ(1) thisQ(2)];
            end
            switch model.teachtype
                case 'mb-teach-mf'
                    thisQH = thisQ + params.replayalpha*(EV-thisQ);
                case 'mf-teach-mb'
                    thisQH = EV + params.replayalpha*(thisQ-EV);
            end
    
            % Determine policy from new hybrid Q value and also what the previous policy would have been
            prob_hybrid = exp([sum(thisQH(1:2)) thisQH(3)]*params.tau) ./ sum(exp([sum(thisQH(1:2)) thisQH(3)]*params.tau));
            prob_mf = exp([sum(thisQ(1:2)) thisQ(3)]*params.tau) ./ sum(exp([sum(thisQ(1:2)) thisQ(3)]*params.tau));
            prob_mb = exp([sum(EV(1:2)) EV(3)]*params.tau) ./ sum(exp([sum(EV(1:2)) EV(3)]*params.tau));
    
            % Compare new value estimates to what value estimates from previous policy (without MB influence)
            switch model.gaintype
                case 'new-old'
                    switch model.teachtype
                        case 'mb-teach-mf'
                            gain = prob_hybrid([1 1 2]).*thisQH - prob_mf([1 1 2]).*thisQH;
                        case 'mf-teach-mb'
                            gain = prob_hybrid([1 1 2]).*thisQH - prob_mb([1 1 2]).*thisQH;
                    end
                case 'old-new'
                    switch model.teachtype
                        case 'mb-teach-mf'
                            gain = prob_mf([1 1 2]).*thisQH - prob_hybrid([1 1 2]).*thisQH;
                        case 'mf-teach-mb'
                            gain = prob_mb([1 1 2]).*thisQH - prob_hybrid([1 1 2]).*thisQH;
                    end
            end

%             oldchoice = double(prob_mf(1)>prob_mf(2));
%             newchoice = double(prob_hybrid(1)>prob_hybrid(2));
%             biggestreplay = find(gain(1:2)==max(gain(1:2)));
%             if length(biggestreplay)==1
%                 if (biggestreplay==1 && d.nV_1(trl)>d.nV_2(trl)) || (biggestreplay==2 && d.nV_1(trl)<d.nV_2(trl))
%                     biggestreplay = 'rewarding';
%                 else
%                     biggestreplay = 'aversive';
%                 end
%                 if oldchoice==newchoice
%                     disp(['TRL ' num2str(trl) ': no policy change (choice = ' num2str(oldchoice) '), biggest replay = ' biggestreplay])
%                 else
%                     disp(['TRL ' num2str(trl) ': change from ' num2str(oldchoice) ' to ' num2str(newchoice) ', biggest replay = ' biggestreplay])
%                 end
%             end
    
            % Make hybrid Q values dependent on which path was replayed
            gainidx = find(gain(1:2) > params.gain);
            thisQH = thisQ;
            switch model.teachtype
                case 'mb-teach-mf' % technically a "model-based" model
                    thisQH(gainidx) = thisQ(gainidx) + params.replayalpha*(EV(gainidx)-thisQ(gainidx));
                case 'mf-teach-mb'
                    thisQH(gainidx) = EV(gainidx) + params.replayalpha*(thisQ(gainidx)-EV(gainidx));
            end
    
            % Recompute policy
            prob = exp([sum(thisQH(1:2)) thisQH(end)]*params.tau) ./ sum(exp([sum(thisQH(1:2)) thisQH(end)]*params.tau)); % collapse across paths
    
%             if length(gainidx)==0
%                 disp(['--- no paths replayed. Choice = ' num2str(double(prob(1)>prob(2)))])
%             elseif length(gainidx)==1
%                 if (d.nV_1(trl)>d.nV_2(trl) && gainidx==1) || (d.nV_1(trl)<d.nV_2(trl) && gainidx==2)
%                     disp(['--- rewarding path replayed. Choice = ' num2str(double(prob(1)>prob(2)))])
%                 else
%                     disp(['--- aversive path replayed. Choice = ' num2str(double(prob(1)>prob(2)))])
%                 end
%             else
%                 disp(['--- BOTH paths replayed. Choice = ' num2str(double(prob(1)>prob(2)))])
%             end


            % Add to table
            d.gain_1(trl) = gain(1);
            d.gain_2(trl) = gain(2);
    
            if d.EV(trl)>1 % should approach
                d.gain_optimal(trl) = gain(1);
                d.gain_suboptimal(trl) = gain(2);
            else % should avoid
                d.gain_optimal(trl) = gain(2);
                d.gain_suboptimal(trl) = gain(1);
            end
    
            d.replayinitiated_1(trl) = gain(1) > params.gain;
            d.replayinitiated_2(trl) = gain(2) > params.gain;

            d.MF_policy(trl) = prob_mf(1);
            d.MB_policy(trl) = prob_mb(1);
            d.MH_policy(trl) = prob_hybrid(1);

            d.Q_1(trl) = thisQ(1);
            d.Q_2(trl) = thisQ(2);
            d.Q_3(trl) = thisQ(3);
    
            d.QH_1(trl) = thisQH(1);
            d.QH_2(trl) = thisQH(2);
            d.QH_3(trl) = thisQH(3);
    
            d.V_1(trl) = EV(1);
            d.V_2(trl) = EV(2);
            d.V_3(trl) = EV(3);
    
            % Make choice
            thisQH = [sum(thisQH(1:2)) thisQH(end)]; % collapse across paths
            prob = exp(thisQH*params.tau) ./ sum(exp(thisQH*params.tau));
    
        end
    
        if all(prob==1) || all(prob==0) || any(isnan(prob))
            warning('Exponentiated number probably too big...? Try making tau smaller')
            errorcatch = true;
        end
    
        if any(prob==0)
            prob(prob==0) = realmin;
            prob(prob==realmin) = 1-realmin;
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
    
        switch MB.calctype
            case {'calculation','experience'}
                if strcmp(MB.calctype,'experience') && strcmp(MB.qtype,'perneg') && ~isempty(V)
                    d.pV_1(trl) = V(nCombos(d.nCombo(trl),1),1);
                    d.pV_2(trl) = V(nCombos(d.nCombo(trl),2),2);
                else
                    if ~isempty(V)
                        d.pV_1(trl) = sum(V(:,1));
                        d.pV_2(trl) = sum(V(:,2));
                    end
                end
            case 'actions'
                switch MB.qtype
                    case 'peraction'
                        d.pV_1(trl) = Q(1);
                        d.pV_2(trl) = Q(2);
                    case 'perprob'
                        d.pV_1(trl) = Q(1,find(pConditions==d.P(trl)));
                        d.pV_2(trl) = Q(2,find(pConditions==d.P(trl)));
                    case 'perneg'
                        d.pV_1(trl) = Q(1,d.nCombo(trl));
                        d.pV_2(trl) = Q(2,d.nCombo(trl));
                    case 'pernegprob'
                        d.pV_1(trl) = Q(1,pConditions==d.P(trl),d.nCombo(trl));
                        d.pV_2(trl) = Q(2,pConditions==d.P(trl),d.nCombo(trl));
                end
        end
    
        % If there was a transition, update perceived value based on outcomes
        if transition>0
    
            outcome = d.Outcome(trl);
    
            % Update values
            switch MB.calctype
                case 'calculation'
                    if any(ismember(MB.pathtype,{'rewarding','aversive'}))
                        V(transition) = V(transition) + params.alpha*(outcome - V(transition));
                    end
                case 'experience'
                    if ~isempty(V)
                        switch MB.qtype
                            case 'perpath'
                                V(transition) = V(transition) + params.alpha*(outcome - V(transition));
                            case 'perneg'
                                V(nCombos(d.nCombo(trl),transition),transition) = V(nCombos(d.nCombo(trl),transition),transition) + params.alpha*(outcome - V(nCombos(d.nCombo(trl),transition),transition));
                        end
                    end
            end
            if transition==1
                statevalues(transition,:) = [d.S1(trl) d.S2(trl) d.S3(trl)];
            elseif transition==2
                statevalues(transition,:) = [d.S4(trl) d.S5(trl) d.S6(trl)];
            end
    
            if strcmp(model.plantype,'hybrid')
                switch MF.calctype
                    case 'experience'
                        switch MF.qtype
                            case 'perpath'
                                Q(transition) = Q(transition) + params.alpha*(outcome - Q(transition));
                            case 'perneg'
                                Q(transition,nCombos(d.nCombo(trl),transition)) = Q(transition,nCombos(d.nCombo(trl),transition)) + params.alpha*(outcome - Q(transition,nCombos(d.nCombo(trl),transition)));
                        end
                end
            end
        end
    
        % Update for either choice, if learning from actions
        switch MF.calctype
            case 'actions'
    
                choice = [];
                if transition>0
                    choice = 1; % approach
                else
                    choice = 2; % avoid
                end
                outcome = d.Outcome(trl);
        
                switch MF.qtype
                    case 'peraction'
                        Q(choice) = Q(choice) + params.alpha*(outcome - Q(choice));
                    case 'perprob'
                        Q(choice,pConditions==d.P(trl)) = Q(choice,pConditions==d.P(trl)) + params.alpha*(outcome - Q(choice,pConditions==d.P(trl)));
                    case 'perneg'
                        Q(choice,d.nCombo(trl)) = Q(choice,d.nCombo(trl)) + params.alpha*(outcome - Q(choice,d.nCombo(trl)));
                    case 'pernegprob'
                        Q(choice,pConditions==d.P(trl),d.nCombo(trl)) = Q(choice,pConditions==d.P(trl),d.nCombo(trl)) + params.alpha*(outcome - Q(choice,pConditions==d.P(trl),d.nCombo(trl)));
                end
            case 'experience'
    
                if ~isempty(Q)
                    outcome = d.Outcome(trl);
                    if transition==0
                        switch MF.qtype
                            case 'perpath'
                                Q(3) = Q(3) + params.alpha*(outcome - Q(3));
                            case 'perneg'
                                Q(3,1) = Q(3,1) + params.alpha*(outcome - Q(3,1));
                        end
                    end
                end
        end
    
        % Add forgetting
        switch model.plantype
            case 'hybrid'
                
                try
                    Q = (1-params.MF_forgetting)*Q; % + (params.MF_forgetting*mean(d.Outcome(1:trl))); 
                catch
                    Q = (1-params.forgetting)*Q;
                end
            
                if ~strcmp(MB.calctype,'calculation') && ~isempty(V)
                    try
                        V = (1-params.MB_forgetting)*V; % + (params.MB_forgetting*mean(d.Outcome(1:trl)));
                    catch
                        V = (1-params.forgetting)*V;
                    end
                else
                    prevstatevalues = (prevstatevalues + statevalues) / trl;
                    try
                        statevalues = (1-params.MB_forgetting)*statevalues + (params.MB_forgetting*prevstatevalues);
                    catch
                        statevalues = (1-params.forgetting)*statevalues + (params.forgetting*prevstatevalues);
                    end
                end
            otherwise
                
                if ~isempty(Q)
                    Q = (1-params.forgetting)*Q;
                end
            
                if ~strcmp(model.calctype,'calculation') && ~isempty(V)
                    V = (1-params.forgetting)*V;
                else
                    prevstatevalues = (prevstatevalues + statevalues) / trl;
                    statevalues = (1-params.forgetting)*statevalues + (params.forgetting*prevstatevalues);
                end
        end
    
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
            T.yhat(trl) = datasample([1 0],1,'replace',false,'weights',pchoice(trl,:));
        end
    end

    err = -nansum(log(T.prob(d.Forced==0)));
    
    % Put into output structure
    output = [];
    output.T = T;
    output.nLL = err;

    % Add counterfactual variables
    if strcmp(model.plantype,'hybrid')

        d.gain_rewarding = nan(size(d,1),1);
        idx = d.nV_1 > d.nV_2; d.gain_rewarding(idx) = d.gain_1(idx);
        idx = d.nV_1 < d.nV_2; d.gain_rewarding(idx) = d.gain_2(idx);
        
        d.gain_aversive = nan(size(d,1),1);
        idx = d.nV_1 < d.nV_2; d.gain_aversive(idx) = d.gain_1(idx);
        idx = d.nV_1 > d.nV_2; d.gain_aversive(idx) = d.gain_2(idx);
        
        d.gain_counterfactual = nan(size(d,1),1);
        d.gain_counterfactual(d.Choice==1) = d.gain_aversive(d.Choice==1);
        d.gain_counterfactual(d.Choice==0) = d.gain_rewarding(d.Choice==0);

        d.replayinitiated_rewarding = nan(size(d,1),1);
        idx = d.nV_1 > d.nV_2; d.replayinitiated_rewarding(idx) = d.replayinitiated_1(idx);
        idx = d.nV_1 < d.nV_2; d.replayinitiated_rewarding(idx) = d.replayinitiated_2(idx);
        
        d.replayinitiated_aversive = nan(size(d,1),1);
        idx = d.nV_1 < d.nV_2; d.replayinitiated_aversive(idx) = d.replayinitiated_1(idx);
        idx = d.nV_1 > d.nV_2; d.replayinitiated_aversive(idx) = d.replayinitiated_2(idx);
        
        d.replayinitiated_counterfactual = nan(size(d,1),1);
        d.replayinitiated_counterfactual(d.Choice==1) = d.replayinitiated_aversive(d.Choice==1);
        d.replayinitiated_counterfactual(d.Choice==0) = d.replayinitiated_rewarding(d.Choice==0);

    end
    
%     disp(['Prediction accuracy = ' num2str(round(mean(d.Choice==T.yhat)*100,2)) '%'])
%     disp(['neg LL              = ' num2str(err)])
%     disp('-------------------------------------------')
end
end