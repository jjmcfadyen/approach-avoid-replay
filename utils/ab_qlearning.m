function [matchrate,greedy_acc,Q,act_greedy] = ab_qlearning(results,alpha)

nTrls = size(results,1);
Q = zeros(2,3,nTrls+1);

Q(1,3,:) = NaN; % assume the Q learner knows the safe value
Q(2,1:2,:) = NaN;
Q(2,3,:) = 1;

act_greedy = nan(nTrls,1);
for trl = 1:nTrls
    
%     % reset q values at the start of blocks
%     if trl > 1 && (results.Block(trl)~=results.Block(trl-1) || results.Practice(trl)~=results.Practice(trl-1))
%         Q(:,:,trl) = zeros(2,3,1);
%         Q(1,3,:) = NaN; % assume the Q learner knows the safe value
%         Q(2,1:2,:) = NaN;
%         Q(2,3,:) = 1;
%     end
    
    % participant choices
    thisOutcome = results.Outcome(trl);
    thisChoice = results.Choice(trl);
    if thisChoice == 1
        thisTransition = results.Transition(trl);
    else
        thisTransition = 3;
    end
    
    % how would a Q learner update values based on this choice & outcome?
    if thisChoice == 1 % can only learn on risky trials
        oldQ = squeeze(Q(:,:,trl));
        newQ = oldQ;
        newQ(thisChoice,thisTransition) = oldQ(thisChoice,thisTransition) + alpha*(thisOutcome - oldQ(thisChoice,thisTransition));
        Q(:,:,trl+1) = newQ;
    else
        Q(:,:,trl+1) = Q(:,:,trl);
    end
    
    % what would the Q learner choose next?
    if trl < nTrls
        
        % take the probability into account
        trlQ = squeeze(Q(:,:,trl+1));
        trlQ(1,1:2) =  trlQ(1,1:2) .* [results.P(trl+1) (1-results.P(trl+1))];
        trlQ_mean = nansum(trlQ,2);
        
        % --- EPSILON-GREEDY
        [~,act_greedy(trl+1)] = max(trlQ_mean);
        
    end
end

%% Get accuracy

greedy_acc = zeros(size(act_greedy,1),1);
greedy_acc(results.EV >= 1 & act_greedy==1) = 1;
greedy_acc(results.EV <= 1 & act_greedy==2) = 1;

% compare to actual behaviour
matchrate = mean(act_greedy == results.Choice);

end