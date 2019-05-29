%% Trials (Practice or Real)

if ispractice
    N = nTrials_practice; % number of weeks
    this_resp_time = practice_resp_time;
    B = 1;
else
    N = nTrials;
    this_resp_time = resp_time;
    B = nBlocks;
end

total_exp_trials = N*nDays;

% Randomly change outcome locations
outcome_locations = [reward_state predator_state]; % col1 = reward, col2 = predator
% if ~ispractice && swap_freq ~= 0
%     average_train_length = swap_freq*total_exp_trials;
%     x_idx = 1:average_train_length*2;
%     train_length_gauss = normpdf(x_idx,average_train_length,swap_sd/2*total_exp_trials);
%     generate = true;
%     while generate
%         this_train = randsample(x_idx,1,true,train_length_gauss);
%         if size(outcome_locations,1) == 1
%             outcome_locations = repmat(outcome_locations,this_train,1);
%         else
%             outcome_locations = [outcome_locations
%                                  repmat(fliplr(outcome_locations(end,1:end)),this_train,1)];
%         end
%         if size(outcome_locations,1) > total_exp_trials
%             outcome_locations = outcome_locations(1:total_exp_trials,:);
%             generate = false;
%         end
%     end
% end

% Compute optimal decisions
this_risk = risk_prob(3);
likely_rewards = [];

% work out probabilistic reward value of each state
state_idx = zeros(length(states),1);
state_idx(states == outcome_locations(1,1)) = reward_mag;
state_idx(states == outcome_locations(1,2)) = -loss_mag;
for i = 1:length(risk_state)
    state_idx(states == risk_state(i)) = this_risk*(-loss_mag);
end
state_idx(4) = (small_reward_prob*small_reward_mag) + ((1-small_reward_prob)*(-starve_mag));

trialinfo =  array2table(nan(total_exp_trials,16), 'VariableNames', ...
    {'block','trial','week','day','predator_prob','forced','move1','state1','rt1','move2','state2','rt2','move3','state3','rt3','health'});

start_state = 1;

baseRectUpper_left = CenterRectOnPointd([0 0 stimulus_size], xCentre-200, yCentre-100);
baseRectUpper_right = CenterRectOnPointd([0 0 stimulus_size], xCentre+200, yCentre-100);
baseRect = CenterRectOnPointd([0 0 stimulus_size], xCentre, yCentre);
baseRect_left = CenterRectOnPointd([0 0 stimulus_size], xCentre-300, yCentre);
top_right_rect = CenterRectOnPointd([0 0 150 150], screenXpixels-100, 125);
top_mid_rect = CenterRectOnPointd([0 0 150 150], xCentre, 125);

t_idx = [this_resp_time-1:-1:0; round((1:this_resp_time)/ifi)]; % for the timer (count every second)

response_coords = [ xCentre-150, yCentre+300
    xCentre,     yCentre+300
    xCentre+150, yCentre+300];

total_rewards = 0;

first_wakefulness_change = [];

trial = 0;
for b = 1:B
    
    if b == 1 && ~ispractice
       
%         DrawFormattedText(window,['We are ready to start the actual experiment.',...
%             '\n\n NOTE! Throughout the experiment, the location of the cheese and the predator',...
%             '\n\n might change on different weeks! Keep track of this.',...
%             '\n\n\n Press ENTER to begin'],'center','center',[0 0 0]);
%         Screen('Flip',window);
%         
%         presstocontinue(enterKey,escapeKey);
        
    end
    
    for week = 1:N
        
        week_reward = 0;
        
        health = datasample([.4 .5 .6],1); % randomly assign percent health to start with
        
        if ispractice
            if week == 1
                predator_wakefulness = 1; % 1 = lowprob, 2 = medprob, 3 = highprob
            elseif week == 2
                predator_wakefulness = 2;
            end
        end
        this_day = 1;
        
        % choose forced choice day(s)
        if forced_choice_freq ~= 0
            if ispractice || (~ispractice && week == 1 && b == 1)
                forced_day_idx = datasample([3:nDays],forced_choice_freq); % pick a random day(s) from the week (but not the first 2 days, for the practice
            else forced_day_idx = datasample(1:nDays,forced_choice_freq);
            end
        else forced_day_idx = 0;
        end

%         % Splash screens
%         DrawFormattedText(window,['Week ' num2str(week) ', Day ' num2str(this_day)],'center','center',[0 0 0]);
%         Screen('Flip',window);
%         WaitSecs(3);
%         
%         if week == 1 && b == 1
%             DrawFormattedText(window,'You will ALWAYS start your exploration from this room:','center',yCentre-200,[0 0 0]);
%             Screen('DrawTexture', window, img(states(start_state)), [], baseRect);
%             Screen('Flip',window);
%             WaitSecs(6);
%         end
%         
%         if week == 1
%             DrawFormattedText(window,['Your starting health is ' num2str(health*100) '%'],'center','center',[0 0 0]);
%         else
%             DrawFormattedText(window,['Your starting health has been reset to ' num2str(health*100) '%'],'center','center',[0 0 0]);
%         end
%         healthbar
%         Screen('Flip',window);
%         WaitSecs(4);
%         
        if ispractice
            if predator_wakefulness == 1
                DrawFormattedText(window,['Today, the predator is ' risk_names{predator_wakefulness} '.',...
                    '\n\n You will only get attacked if you walk straight into it.'],'center','center',[0 0 0]);
            elseif predator_wakefulness == 2
                DrawFormattedText(window,['Today, the predator is ' risk_names{predator_wakefulness} '.',...
                    '\n\n If you enter it’s path, it might attack you!'],'center','center',[0 0 0]);
            elseif predator_wakefulness == 3
                DrawFormattedText(window,['Today, the predator is ' risk_names{predator_wakefulness} '.',...
                    '\n\n Be careful! It is likely to attack if you enter it’s path!'],'center','center',[0 0 0]);
            end
            Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
            Screen('Flip',window);
            WaitSecs(4);
        elseif ~ispractice && week == 1 && b == 1
            predator_wakefulness = randi(length(risk_prob)); % randomly pick the predator probability for this day
%             DrawFormattedText(window,'Remember to always check the top right icon to see how awake the predator is.','center','center',[0 0 0]);
%             Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
%             Screen('Flip',window);
%             WaitSecs(4);
        end
%         
%         if b == 1 && week == 1
%             DrawFormattedText(window,'Let’s start!','center','center',[0 0 0]);
%             healthbar;
%             Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
%             Screen('Flip',window);
%             WaitSecs(2);
%         end
        

        
        for this_day = 1:nDays

            if ispractice
                if week == 2 && this_day > 1
                    old_wakefulness = predator_wakefulness;
                    new_wakefulness = datasample([2 3],1); % randomly pick whether the predator is fully awake or half awake
                    if new_wakefulness ~= old_wakefulness
                        first_wakefulness_change = [first_wakefulness_change, 1];
                    end
                    predator_wakefulness = new_wakefulness;
                end
            else
                if b == 1 && week == 1 && this_day == 1
                else predator_wakefulness = randi(length(risk_prob)); % randomly pick the predator probability for this day
                end
            end
            
            if health > 0
                
                if this_day > 1
                    DrawFormattedText(window,['Week ' num2str(week) ', Day ' num2str(this_day)],'center','center',[0 0 0]);
                    Screen('Flip',window);
                    WaitSecs(3);
                end
                
                if length(first_wakefulness_change) == 1
%                     DrawFormattedText(window,['ATTENTION! The predator’s wakefulness has changed from ' risk_names{old_wakefulness},...
%                         ' to ' risk_names{new_wakefulness} '. \n\n See the icon in the top right has changed.'],'center','center',[0 0 0]);
%                     Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
%                     Screen('Flip',window);
%                     WaitSecs(4);
                    first_wakefulness_change = false;
                elseif length(first_wakefulness_change) == 2
%                     DrawFormattedText(window,'Remember to always check the top right icon to see how awake the predator is.','center','center',[0 0 0]);
%                     Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
%                     Screen('Flip',window);
%                     WaitSecs(4);
                end
                
                if any(this_day == forced_day_idx)
                    if ispractice || this_day == 1 && week == 1 && b == 1
                        DrawFormattedText(window,...
                            '** On this day, you must go to a SPECIFIC room **',...
                            '\n\n Enter in the moves that will get you there (despite the predator/cheese)',...
                            '\n\n There is no time limit.','center','center',[0 0 0]);
                        Screen('Flip',window);
                        WaitSecs(4);
                    else
                        DrawFormattedText(window,...
                            '** On this day, you must go to a SPECIFIC room **','center','center',[0 0 0]);
                        Screen('Flip',window);
                        WaitSecs(3);
                    end
                    
                    % Choose the least explored path
                    visited_states = [trialinfo.state1, trialinfo.state2, trialinfo.state3];
                    visited_states = visited_states(~isnan(visited_states(:,1)),:);
                    state_freq = [];
                    for st = 1:length(states)
                        state_freq(st,:) = [states(st), sum(visited_states(:) == states(st))];
                    end
                    state_freq(state_freq(:,1) == start_state,:) = []; % ignore starting state
                    min_state_freq = state_freq(state_freq(:,2) == min(state_freq(:,2))); % find least-visited states
                    min_paths = zeros(size(all_paths,1),size(all_paths,2));
                    for p = 1:size(all_paths,1)
                        for st = 1:length(min_state_freq)
                            min_paths(p,all_paths(p,:) == min_state_freq(st)) = 1;
                        end
                    end
                    least_explored_path = find(sum(min_paths,2) == max(sum(min_paths,2)));
                    if length(least_explored_path) > 1
                        least_explored_path = least_explored_path(randi(length(least_explored_path))); % randomly pick
                    end
                    least_explored_path = all_paths(least_explored_path,:);
                    
                end
                
                if health > 1
                    health = 1; % can't be greater than 100%
                end
                
                trial = trial + 1;
                trialinfo.trial(trial) = trial;
                trialinfo.week(trial) = week;
                trialinfo.day(trial) = this_day;
                
                if any(this_day == forced_day_idx)
                    trialinfo.forced(trial) = 1;
                else trialinfo.forced(trial) = 0;
                end
                
                if ~ispractice && swap_freq ~= 0
                    predator_state = outcome_locations(trial,2);
                    reward_state = outcome_locations(trial,1);
                end
                
                % display trial information
                disp(['Block ' num2str(b) ', trial ' num2str(trial) ': predator is ' risk_names{predator_wakefulness}])
                
                answer_complete = false;
                finish_q = false;
                moves_entered = [0 0 0];
                timer = this_resp_time;
                timeout = false;
                f = 1; % frame rate counter for timer
                keyCode = zeros(1,256);
                while ~answer_complete
                    
                    % set keyboard parameters (look for a CHANGE in whether it was
                    % pressed/released - otherwise it think a held key is multiple
                    % responses)
                    old_keyCode = keyCode;
                    
                    % see what move we're up to
                    nMove = find(moves_entered == 0);
                    if length(nMove) > 1
                        nMove = nMove(1);
                    end
                    
                    % Draw the question again every frame
                    Screen('DrawTexture', window, img(states(start_state)), [], baseRectUpper_left);
                    DrawFormattedText(window,'to','center',yCentre-100,[0 0 0]);
                    if any(this_day == forced_day_idx)
                        Screen('DrawTexture', window, img(states(least_explored_path(end))), [], baseRectUpper_right);
                    else
                        DrawFormattedText(window,'?',xCentre+200,yCentre-100,[0 0 0]);
                    end
                    Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
                    healthbar; % runs external script
                    
                    % Draw timer (but not for forced-choice trials)
                    if ~any(this_day == forced_day_idx)
                        if any(t_idx(2,:) == f)
                            timer = t_idx(1,t_idx(2,:) == f); % for the timer (count every second)
                            if timer == 0
                                timeout = true;
                                answer_complete = true;
                            end
                        end
                        DrawFormattedText(window,['[ ' num2str(timer) ' seconds remaining ]'],'center',yCentre+150,[0 0 1]);
                    end
                    
                    [~, ~, keyCode, ~] = KbCheck;
                    keyDiff = keyCode - old_keyCode; % 1 = newly pressed, -1 = newly released
                    keyDiff(keyDiff == -1) = 0; % ignore the 'newly released' keys
                    if keyDiff(escapeKey)
                        sca
                        return
                    elseif ~isempty(nMove) && (keyDiff(upKey) || keyDiff(downKey) || keyDiff(leftKey) || keyDiff(rightKey))
                        this_move = find(key_map == find(keyDiff));
                        moves_entered(nMove) = this_move;
                        %                         disp(['Move ' num2str(nMove) ' = ' action_labels{this_move}])
                        %                         disp([num2str(response_coords(nMove,1)) ', ' num2str(response_coords(nMove,2))])
                        if sum(moves_entered == 0) == 0
                            finish_q = true;
                            DrawFormattedText(window,'Answer submitted!','center',yCentre+350,[0 1 0]);
                        end
                        rt = GetSecs-t0;
                        
                        % update table
                        if nMove == 1
                            trialinfo.move1(trial) = this_move;
                            trialinfo.rt1(trial) = rt;
                        elseif nMove == 2
                            trialinfo.move2(trial) = this_move;
                            trialinfo.rt2(trial) = rt;
                        elseif nMove == 3
                            trialinfo.move3(trial) = this_move;
                            trialinfo.rt3(trial) = rt;
                        end
                    end
                    
                    % Put in the response placeholders
                    if ~finish_q
                        arrow_colour = [0 0 0];
                    else arrow_colour = [0 1 0];
                    end
                    DrawFormattedText(window,'Put in your 3 moves','center',yCentre+200,[0 0 0]);
                    for r = 1:3
                        if moves_entered(r) == 0
                            DrawFormattedText(window,'-',response_coords(r,1),response_coords(r,2),arrow_colour); % placeholders for arrow responses
                        else
                            DrawFormattedText(window,response_text{moves_entered(r)},response_coords(r,1),response_coords(r,2),arrow_colour); % arrow response
                        end
                    end
                    
                    vbl = Screen('Flip',window);
                    if f == 1
                        t0 = vbl;
                    end
                    f = f + 1;
                    
                    if finish_q == true
                        WaitSecs(2);
                        answer_complete = true;
                    end
                    
                end
                
                if timeout && ~any(this_day == forced_day_idx)
                    
                    healthbar
                    DrawFormattedText(window,'You ran out of time!','center','center',[1 0 0]);
                    healthbar;
                    Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
                    Screen('Flip',window);
                    WaitSecs(3);

                    DrawFormattedText(window,['You went hungry and lost ' num2str(starve_mag*100) '% health'],'center','center',[1 0 0]);
                    health = health - starve_mag;
                    if health < 0
                        health = 0;
                    end
                    healthbar;
                    Screen('DrawTexture', window, threatprob_stimuli(predator_wakefulness), [], top_right_rect);
                    Screen('Flip',window);
                    WaitSecs(3);
                    
                    % check if died
                    if health == 0
                        DrawFormattedText(window,'You died of starvation','center','center',[1 0 0]);
                        Screen('Flip',window);
                        WaitSecs(3);
                    end
                    
                else
                    
                    % Give feedback
                    DrawFormattedText(window,'Let´s follow your chosen path...','center',yCentre-200,[0 0 0]);
                    Screen('DrawTexture', window, img(states(start_state)), [], baseRect); % draw start state for question
                    Screen('Flip',window);
                    WaitSecs(2);
                    states_transitioned = [0 0 0];
                    calculate_risk = 0;
                    for r = 1:3 % show the transitions for each move...
                        
                        if health > 0
                            
                            this_action = moves_entered(r);
                            
                            % get start and end state based on the moves they entered
                            if r == 1
                                state1 = start_state;
                            end
                            
                            [~, state2] = state_transition_withhealth(this_action,state1,transitions,health,this_day,params,1);
                            
                            if any(this_day == forced_day_idx)
                                if state2 ~= least_explored_path(r) % if they went the wrong way
                                    
                                    WaitSecs(.5);
                                    healthbar;
                                    DrawFormattedText(window,'You went the wrong way...','center',yCentre-200,[1 0 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('Flip',window);
                                    WaitSecs(2);
                                    
                                    % move them back to previous state
                                    if state1 ~= state2
                                        reverse_action = transitions(transitions(:,1) == state2 & transitions(:,3) == state1,2);
                                        [~, state2] = state_transition_withhealth(reverse_action,state1,transitions,health,this_day,params,0);
                                    end
                                    
                                    % correct action
                                    correct_action = transitions(transitions(:,1) == state1 & transitions(:,3) == least_explored_path(r),2);
                                    [~, state2] = state_transition_withhealth(correct_action,state1,transitions,health,this_day,params,1);
                                    
                                end
                            end
                            
                            states_transitioned(r) = state2;
                            if r == 1
                                trialinfo.state1(trial) = state2;
                            elseif r == 2
                                trialinfo.state2(trial) = state2;
                            elseif r == 3
                                trialinfo.state3(trial) = state2;
                            end
                            
                            if any(state2 == risk_state)
                                calculate_risk = randsample([0 1],1,true,[1-risk_prob(predator_wakefulness) risk_prob(predator_wakefulness)]);
                                disp(['Risky state = ' num2str(calculate_risk)])
                            elseif state2 == predator_state
                                calculate_risk = 1;
                            else calculate_risk = 0;
                            end
                            
                            if calculate_risk == 1
                                if state2 == predator_state || (any(state2 == risk_state) && ~any(this_day == forced_day_idx))
                                  
                                    % only have predators attack on risk
                                    % states if it's NOT a forced choice
                                    % trial (the point of these is to get people to the end
                                    % of the path, so we don't want them to
                                    % die on the way)
                                    
                                    % calculate loss
                                    if loss_mag == Inf
                                        total_loss = health;
                                        loss_string = 'all of your';
                                    else
                                        total_loss = loss_mag;
                                        loss_string = [num2str(total_loss*100) '%'];
                                    end
                                    % draw attack screen with health BEFORE attack
                                    WaitSecs(.5);
                                    healthbar;
                                    DrawFormattedText(window,'You were attacked!','center',yCentre-200,[1 0 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('DrawTexture', window, outcome_stimuli(2), [], baseRect_left);
                                    DrawFormattedText(window,['You lost ' loss_string ' health'],'center',yCentre+200,[1 0 0]);
                                    Screen('Flip',window);
                                    WaitSecs(2);
                                    % draw attack screen with health AFTER attack
                                    health = health-total_loss;
                                    if health < 0
                                        health = 0;
                                    end
                                    healthbar;
                                    DrawFormattedText(window,'You were attacked!','center',yCentre-200,[1 0 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('DrawTexture', window, outcome_stimuli(2), [], baseRect_left);
                                    DrawFormattedText(window,['You lost ' loss_string ' health'],'center',yCentre+200,[1 0 0]);
                                    Screen('Flip',window);
                                    WaitSecs(1.5);
                                    % check if died
                                    if health == 0
                                        DrawFormattedText(window,'You were killed by the predator','center','center',[1 0 0]);
                                        Screen('Flip',window);
                                        WaitSecs(3);
                                    end
                                end
                            elseif state2 == reward_state
                                WaitSecs(.5);
                                % draw health BEFORE reward
                                healthbar;
                                if ~any(this_day == forced_day_idx)
                                    DrawFormattedText(window,'You found cheese!','center',yCentre-200,[0 1 0]);
                                    DrawFormattedText(window,['You gained ' num2str(reward_mag*100) '% health'],'center',yCentre+200,[0 1 0]);
                                else
                                    DrawFormattedText(window,'You would have found cheese!','center',yCentre-200,[0 1 0]);
                                    DrawFormattedText(window,['You would have gained ' num2str(reward_mag*100) '% health'],'center',yCentre+200,[0 1 0]);
                                end
                                Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                Screen('DrawTexture', window, outcome_stimuli(1), [], baseRect_left);
                                Screen('Flip',window);
                                WaitSecs(2);
                                if ~any(this_day == forced_day_idx)
                                    % draw health AFTER reward
                                    health = health + reward_mag;
                                    if health > 1
                                        health = 1;
                                    end
                                    healthbar;
                                    DrawFormattedText(window,'You found cheese!','center',yCentre-200,[0 1 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('DrawTexture', window, outcome_stimuli(1), [], baseRect_left);
                                    DrawFormattedText(window,['You gained ' num2str(reward_mag*100) '% health'],'center',yCentre+200,[0 1 0]);
                                    Screen('Flip',window);
                                    WaitSecs(2);
                                end
                            end
                            
                            if r == 3 && state2 ~= predator_state && ~any(risk_state == state2) && state2 ~= reward_state
                                % crumbs or not?
                                found_crumbs = 0;
                                if small_reward_mag > 0 && small_reward_prob > 0
                                    found_crumbs = randsample([0 1],1,true,[1-small_reward_prob small_reward_prob]);
                                end
                                
                                WaitSecs(0.5);
                                if found_crumbs == 0
                                    % draw feedback BEFORE health loss
                                    DrawFormattedText(window,['There is nothing here.',...
                                        '\n\n You went hungry and lost ' num2str(starve_mag*100) '% health'],'center',yCentre+200,[1 0 0]);
                                    healthbar;
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('Flip',window);
                                    WaitSecs(3);
                                    % draw feedback AFTER health loss
                                    DrawFormattedText(window,['There is nothing here.',...
                                        '\n\n You went hungry and lost ' num2str(starve_mag*100) '% health'],'center',yCentre+200,[1 0 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    health = health - starve_mag;
                                    healthbar;
                                    Screen('Flip',window);
                                    WaitSecs(3);
                                    % check if died
                                    if health == 0
                                        DrawFormattedText(window,'You died of starvation','center','center',[1 0 0]);
                                        Screen('Flip',window);
                                        WaitSecs(3);
                                    end
                                elseif found_crumbs == 1
                                    % draw feedback BEFORE health loss
                                    DrawFormattedText(window,['There wasn´t much down this path...',...
                                        '\n\n But luckily you found some crumbs and gained ' num2str(small_reward_mag*100) '% health'],'center',yCentre+200,[0 0 0]);
                                    healthbar;
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('DrawTexture', window, outcome_stimuli(3), [], baseRect_left);
                                    Screen('Flip',window);
                                    WaitSecs(3);
                                    % draw feedback AFTER health loss
                                    DrawFormattedText(window,['There wasn´t much down this path...',...
                                        '\n\n But luckily you found some crumbs and gained ' num2str(small_reward_mag*100) '% health'],'center',yCentre+200,[1 0 0]);
                                    Screen('DrawTexture', window, img(states(state2)), [], baseRect);
                                    Screen('DrawTexture', window, outcome_stimuli(3), [], baseRect_left);
                                    health = health - small_reward_mag;
                                    healthbar;
                                    Screen('Flip',window);
                                    WaitSecs(3);
                                end
                            end
                            
                            state1 = state2; % update start state for next move
                        end
                        
                    end
                    
                    % update table
                    trialinfo.health(trial) = health;
                    trialinfo.state1(trial) = states_transitioned(1);
                    trialinfo.state2(trial) = states_transitioned(2);
                    trialinfo.state3(trial) = states_transitioned(3);
                    
                end
                
                if health == 0
                    days_survived = this_day-1;
                else
                    days_survived = nDays;
                end
                
            end
        end
        
        if health == 0
            DrawFormattedText(window,['End of week ' num2str(week),...
                '\n\n You managed to survive for ' num2str(days_survived) ' days and then died on day ' num2str(days_survived+1) '.',...
                '\n\n\n Money earnt this week: £0',...
                '\n\n Money earnt so far: £' num2str(total_rewards)],'center','center',[0 0 0]);
            Screen('Flip',window);
            WaitSecs(4);
        elseif health > 0 && health < 1
            total_rewards = total_rewards + 1;
            week_reward = week_reward + 1;
            DrawFormattedText(window,['End of week ' num2str(week),...
                '\n\n You managed to survive for ' num2str(days_survived) ' days and ended with ' num2str(health*100) '% health.',...
                '\n\n\n Money earnt this week: £' num2str(week_reward),...
                '\n\n Money earnt so far: £' num2str(total_rewards)],'center','center',[0 0 0]);
            Screen('Flip',window);
            WaitSecs(4);
        elseif health == 1
            total_rewards = total_rewards + 2;
            week_reward = week_reward + 2;
            DrawFormattedText(window,['End of week ' num2str(week),...
                '\n\n You managed to survive for ' num2str(days_survived) ' days and maintained 100% health!',...
                '\n\n\n Money earnt this week: £' num2str(week_reward),...
                '\n\n Money earnt so far: £' num2str(total_rewards)],'center','center',[0 0 0]);
            Screen('Flip',window);
            WaitSecs(4);
        end
        
        if week < N
            DrawFormattedText(window,'Respawning for new week...','center','center',[0 0 0]);
            Screen('Flip',window);
            WaitSecs(3);
        end
        
    end
    
    if ispractice
        practice_trialinfo = trialinfo;
    end
    
    if ispractice
        DrawFormattedText(window,['End of practice!',...
            '\n\n You earnt £' num2str(total_rewards) ' in total.',...
            '\n\n\n Press ENTER to start the real experiment.'],...
            'center','center',[0 0 0]);
    else
        DrawFormattedText(window,['End of block ' num2str(b),...
            '\n\n You have earnt £' num2str(total_rewards) ' in total.',...
            '\n\n\n Press ENTER to start the next block.'],...
            'center','center',[0 0 0]);
    end
    
    Screen('Flip',window);
    presstocontinue(enterKey,escapeKey);
    
end

if ~ispractice
    DrawFormattedText(window,'End of experiment! Thank you for participating','center','center',[0 0 0]);
    Screen('Flip',window)
    WaitSecs(4);
    sca;
    return
end
