function [err, output, d] = sim_model(d,model,vec)
% [err, model] = sim_model(d,model,vec)
% 
% d = subject data table
% model = model specification (single model from output of set_models --> e.g., model = models(1))
% vec = vector of parameter values (in order specified in 'model' structure)

randtype = 'both'; % 'one' or 'both'

%% Check input

nTrls = size(d,1);

if nargin==2
    params = [];
else
    vec = vec(:)'; % make it a row
end

if ~isfield(model.params,'alpha') % learning rate for which path is rewarding and which is aversive
    params.alpha = 1;
else
    params.alpha = model.params.alpha.val;
end
if ~isfield(model,'type') % more details about which path is computed (e.g., 1, 2, rewarding, aversive)
    params.type = '';
else
    params.type = model.type;
end
if ~isfield(model.params,'threshold') % threshold at which you'll approach (i.e., the subjective value of avoiding)
    params.threshold = 1;
else
    params.threshold = model.params.threshold.val;
end
if ~isfield(model.params,'rewthreshold') % threshold at which you'll approach (i.e., the subjective value of avoiding)
    params.rewthreshold = 0;
else
    params.rewthreshold = model.params.rewthreshold.val;
end
if ~isfield(model.params,'lossthreshold') % threshold at which you'll approach (i.e., the subjective value of avoiding)
    params.lossthreshold = 0;
else
    params.lossthreshold = model.params.lossthreshold.val;
end
if ~isfield(model.params,'tau') % inverse temperature parameter for softmax
    params.tau = 1;
else
    params.tau = model.params.tau.val;
end
if ~isfield(model,'nThresh') % whehter you compare all path values to a single number, or if you have different numbers depending on whether a path is positive or negative
    params.nThresh = 1;
else
    params.nThresh = model.nThresh;
end
if ~isfield(model,'qlearner')
    qlearner = false;
elseif isempty(model.qlearner)
    qlearner = false;
else
    qlearner = model.qlearner;
end
if ~isfield(model,'neglearn')
    neglearn = false;
elseif isempty(model.neglearn)
    neglearn = false;
else
    neglearn = model.neglearn;
end

if nargin==3 % parameters are being entered in optimisation function
    
    pnames = fieldnames(model.params);
    np = length(pnames);
    freeidx = zeros(1,length(pnames));
    for p = 1:np
        freeidx(p) = isnan(model.params.(pnames{p}).val);
    end
    
    pnames = pnames(freeidx==1);
    np = length(pnames);

    for p = 1:np
        params.(pnames{p}) = vec(p);
    end
end

% Check values are within specified ranges
if params.alpha < 0 || params.alpha > 1 || params.tau <= 0
    error('A parameter is out of bounds')
end

% Check if choice data exists or not
existingdata = false;
colnames = d.Properties.VariableNames;
if ~any(contains(colnames,'Choice')) || ~any(contains(colnames,'Transition')) || ~any(contains(colnames,'Outcome'))
    d.Choice = nan(size(d,1),1);
    d.Transition = nan(size(d,1),1);
    d.Outcome = nan(size(d,1),1);
else
    existingdata = true;
end

%% Estimate choices using model parameters

if ~neglearn
    old_nV = [0 0]; % known value of each path
else
    old_nV = zeros(3,2); % rows = where the odd state is, columns = path 1 and path 2
    nCombos = fliplr(combvec(1:3,1:3)');
end

pchoice = nan(nTrls,2); % cols: approach, avoid
errorcatch = false;
for trl = 1:nTrls

    % Get position of odd rule ("neg states"), if applicable
    if neglearn
        negpos = nCombos(d.nCombo(trl),:);
    end

    % If forced-choice, set prediction to NaN
    if d.Forced(trl)>0
        pchoice(trl,1) = NaN;
    else

        % Determine value
        switch model.fun
            case 'twopath'
                if qlearner
                    if neglearn
                        thisnv = [old_nV(negpos(1),1) old_nV(negpos(2),2)];
                        V = [sum(thisnv.*[d.P(trl) 1-d.P(trl)]) params.threshold];
                    else
                        V = [sum(old_nV.*[d.P(trl) 1-d.P(trl)]) params.threshold];
                    end
                else
                    V = [d.EV(trl) params.threshold];
                end
            case 'random'
                V = [params.threshold 0]; % threshold indicates general preference for approach (more positive) vs avoid (more negative)
            case 'onepath'

                % which path's value to consider
                path = [];
                if strcmp(params.type,'random')
                    path = [1 2];
                elseif contains(params.type,'1')
                    path = 1;
                elseif contains(params.type,'2')    
                    path = 2;
                elseif strcmp(params.type,'rewarding') || strcmp(params.type,'aversive')
                    switch params.type
                        case 'rewarding'
                            path = find(sum(old_nV,1)==max(sum(old_nV,1)));
                        case 'aversive'
                            path = find(sum(old_nV,1)==min(sum(old_nV,1)));
                    end
                    if length(path)>1
                        if strcmp(randtype,'one')
                            path = randi(2);
                        else
                            path = [1 2];
                        end
                    end
                else
                    error('Could not determine which path is to be calculated')
                end

                V = nan(length(path),2);
                for p = 1:length(path)

                    P = d.P(trl);
                    if path(p)==2
                        P = 1-P;
                    end
    
                    if qlearner
                        if ~neglearn
                            EV = old_nV(path(p))*P;
                        else
                            thisnv = [old_nV(negpos(1),1) old_nV(negpos(2),2)];
                            EV = thisnv(path(p))*P;
                        end
                    else
                        EV = d.(['nV_' num2str(path(p))])(trl)*P;
                    end
    
                    if params.nThresh==1
                        V(p,:) = [EV params.threshold];
                    elseif params.nThresh==2
                        if EV>0
                            V(p,:) = [EV params.rewthreshold]; % threshold used for rewarding paths
                        else
                            V(p,:) = [EV params.lossthreshold]; % threshold used for aversive paths
                        end
                    end
                end
        end

        % Determine choice probabilities using softmax
        prob = nan(size(V,1),2);
        for p = 1:size(V,1)
            prob(p,:) = exp(V(p,:)*params.tau) ./ sum(exp(V(p,:)*params.tau));
        end
        if size(prob,1)>1
            switch randtype
                case 'both'
                    pchoice(trl,:) = mean(prob);
                case 'one'
                    pchoice(trl,:) = prob(randi(2),:);
            end
        else
            pchoice(trl,:) = prob;
        end

        if all(pchoice(trl,:)==1) || all(pchoice(trl,:)==0) || any(isnan(pchoice(trl,:)))
%             warning('Exponentiated number probably too big...? Try making tau smaller')
            errorcatch = true;
        end

        if any(pchoice(trl,:)==0)
            pchoice(trl,pchoice(trl,:)==0) = realmin;
            pchoice(trl,pchoice(trl,:)==realmin) = 1-realmin;
        end

    end

    % Update knowledge of path values
    if existingdata
        transition = d.Transition(trl);
    else
        if d.Forced(trl)==0
            thischoice = datasample([1 0],1,'replace',false,'weights',pchoice(trl,:)); % choose action with some randomness
        else
            thischoice = 1;
        end
        if thischoice==1
            transition = datasample([1 2],1,'weights',[d.P(trl) 1-d.P(trl)]);
            d.Outcome(trl) = d.(['nV_' num2str(transition)])(trl);
        else
            transition = 0;
            d.Outcome(trl) = 1;
        end
        d.Transition(trl) = transition;
        d.Choice(trl) = thischoice;
    end
    if transition>0
        if ~neglearn
            old_nV(transition) = old_nV(transition) + params.alpha*(d.Outcome(trl) - old_nV(transition));
        else
            old_nV(negpos(transition),transition) = old_nV(negpos(transition),transition) + params.alpha*(d.Outcome(trl) - old_nV(negpos(transition),transition));
        end
    end
end

%% Create output table (averaged over several iterations)

if errorcatch

    err = Inf;
    output = [];

else

    nIterations = 1; % to account for randomness
    
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
    
    err = -nansum(log(T.prob(idx)));
    
    % Put into output structure
    output = [];
    output.T = T;
    output.nLL = err;
    
    % disp(['Prediction accuracy = ' num2str(round(mean(d.Choice==yhat)*100,2)) '%'])
    % disp(['neg LL              = ' num2str(err)])
    % disp(['EV threshold        = ' num2str(params.threshold)])
    % disp(['Learning rate       = ' num2str(params.alpha)])
    % disp('-------------------------------------------')
end
end