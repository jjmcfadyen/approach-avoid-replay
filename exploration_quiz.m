% Instruction screen
DrawFormattedText(window,['Now that you have explored the map, we will do a short quiz.',...
                          '\n\n\n We will show you a START room and a FINISH room.',...
                          '\n\n You will then enter in the series of moves (UP, DOWN, LEFT, RIGHT) that would get you from the start to the finish',...
                          '\n\n\n Press ENTER to begin'],...
                          'center',yCentre);
Screen('Flip',window);
presstocontinue(enterKey,escapeKey);

nQ = size(test_series,1);
min_accuracy = .95;

test_series = test_series(datasample(1:nQ,nQ,'replace',false),:); % randomly shuffle question order

% Question display parameters
start_pos = screenXpixels * .25;
finish_pos = screenXpixels * .75;

% Ask questions
test_answers =  array2table(nan(size(test_series,1),9), 'VariableNames', ...
{'start_state','finish_state','response_move1','response_move2','response_move3','answer_move1','answer_move2','answer_move3','accuracy'});
for q = 1:nQuizTrials

    start_state = test_series(q,1);
    finish_state = test_series(q,2);

    start_rect = CenterRectOnPointd([0 0 300 300], start_pos, yCentre-100);
    finish_rect = CenterRectOnPointd([0 0 300 300], finish_pos, yCentre-100);

    % run
    response_coords = [xCentre-150, yCentre+300
                       xCentre,     yCentre+300
                       xCentre+150, yCentre+300];

    % Draw the question (start & finish stimuli)
    DrawFormattedText(window,'If you started here...',start_pos-120,yCentre-300,[0 0 0]); % try doing this with centering within rectangles instead?
    Screen('DrawTexture', window, img(states(start_state)), [], start_rect);

    DrawFormattedText(window,'... how would you get here?',finish_pos-160,yCentre-300,[0 0 0]);
    Screen('DrawTexture', window, img(states(finish_state)), [], finish_rect);

    DrawFormattedText(window,'Put in your 3 moves. To submit, press ENTER. To clear, press C.','center',yCentre+200,[0 0 0]);

    % Put in the response placeholders
    for r = 1:3
        DrawFormattedText(window,'-',response_coords(r,1),response_coords(r,2)); % placeholders for arrow responses
    end

    Screen('Flip',window);

    answer_complete = false;
    finish_q = false;
    moves_entered = [0 0 0];
    while ~answer_complete

        % see what move we're up to
        nMove = find(moves_entered == 0);
        if length(nMove) > 1
            nMove = nMove(1);
        end

        % Draw the question (start & finish stimuli) again every frame
        DrawFormattedText(window,'If you started here...',start_pos-120,yCentre-300,[0 0 0]); % try doing this with centering within rectangles instead?
        Screen('DrawTexture', window, img(states(start_state)), [], start_rect);

        DrawFormattedText(window,'... how would you get here?',finish_pos-160,yCentre-300,[0 0 0]);
        Screen('DrawTexture', window, img(states(finish_state)), [], finish_rect);

        DrawFormattedText(window,'Put in your 3 moves. To submit, press ENTER. To clear, press C.','center',yCentre+200,[0 0 0]);

        [secs, keyCode, deltaSecs] = KbPressWait;
        if keyCode(escapeKey)
            sca
            return
        elseif ~isempty(nMove) && (keyCode(upKey) || keyCode(downKey) || keyCode(leftKey) || keyCode(rightKey))
            this_move = find(key_map == find(keyCode));
            moves_entered(nMove) = this_move;
            disp(['Move ' num2str(nMove) ' = ' action_labels{this_move}])
            disp([num2str(response_coords(nMove,1)) ', ' num2str(response_coords(nMove,2))])
        elseif keyCode(enterKey) && sum(moves_entered == 0) == 0
            finish_q = true;
            DrawFormattedText(window,'Answer submitted!','center',yCentre+350,[0 1 0]);
        elseif keyCode(enterKey) && any(moves_entered == 0)
            DrawFormattedText(window,'You must enter 3 moves before you can submit!','center',yCentre+400,[1 0 0]);
        elseif keyCode(clearKey)
            moves_entered = [0 0 0];
        end

        % Put in the response placeholders
        if ~finish_q
            arrow_colour = [0 0 0];
        else arrow_colour = [0 1 0];
        end
        for r = 1:3
            if moves_entered(r) == 0
                DrawFormattedText(window,'-',response_coords(r,1),response_coords(r,2),arrow_colour); % placeholders for arrow responses
            else
                DrawFormattedText(window,response_text{moves_entered(r)},response_coords(r,1),response_coords(r,2),arrow_colour); % arrow response
            end
        end

        Screen('Flip',window);

        if finish_q == true
            WaitSecs(2);
            answer_complete = true;
        end

    end

    % fill in table
    test_answers.start_state(q) = test_series(q,1);
    test_answers.finish_state(q) = test_series(q,2);
    test_answers.answer_move1(q) = test_series(q,3);
    test_answers.answer_move2(q) = test_series(q,4);
    test_answers.answer_move3(q) = test_series(q,5);
    test_answers.response_move1(q) = moves_entered(1);
    test_answers.response_move2(q) = moves_entered(2);
    test_answers.response_move3(q) = moves_entered(3);
    test_answers.accuracy(q) = mean(test_series(1,3:end) == moves_entered); % percentage of the moves correctly entered

    % Preview their chosen path
    baseRect = CenterRectOnPointd([0 0 300 300], xCentre, yCentre);

    DrawFormattedText(window,'Let´s follow your chosen path...','center',yCentre-200,[0 0 0]);
    Screen('DrawTexture', window, img(states(test_series(q,1))), [], baseRect); % draw start state for question
    Screen('Flip',window);
    WaitSecs(2);
    states_transitioned = [0 0 0];
    for r = 1:3 % show the transitions for each move...

        this_action = moves_entered(r);

        % get start and end state based on the moves they entered
        if r == 1
            state1 = test_series(q,1);
        end

        [~, state2] = state_transition(this_action,state1,transitions,params,1);

        states_transitioned(r) = state2;
        state1 = state2; % update start state for next move

    end

    % Give correct/incorrect feedback
    if state1 == test_series(q,2) % they ended up in the correct state
        DrawFormattedText(window,'Great work! You found your way!','center',yCentre-200,[0 1 0]);
        Screen('DrawTexture', window, img(states(state1)), [], baseRect);
        Screen('Flip',window)
        WaitSecs(2)
    else
        DrawFormattedText(window,'Oh no, it looks like you lost your way...','center',yCentre-200,[1 0 0]);
        Screen('DrawTexture', window, img(states(state1)), [], baseRect);
        Screen('Flip',window)
        WaitSecs(2)

        % show the correct path
        DrawFormattedText(window,'This was the correct path...','center',yCentre-200,[0 0 0]);
        Screen('DrawTexture', window, img(states(test_series(q,1))), [], baseRect); % draw start state for question
        Screen('Flip',window);
        WaitSecs(2);

        for r = 1:3 % show the CORRECT transitions for each move...

            this_action = test_series(q,2+r);

            % get start and end state based on the moves they entered
            if r == 1
                state1 = test_series(q,1);
            end

            [~, state2] = state_transition(this_action,state1,transitions,params,1);

            states_transitioned(r) = state2;
            state1 = state2; % update start state for next move

        end

    end

end

% Give quiz result
quiz_accuracy = round(nanmean(test_answers.accuracy*100),2);
if quiz_accuracy/100 < quiz_threshold
    DrawFormattedText(window,['You scored ' num2str(quiz_accuracy) '% on the quiz',...
        '\n\n We need you to score at least ' num2str(round(quiz_threshold*100)) '%',...
        '\n\n Let´s go back and explore some more.',...
        '\n\n\n Press ENTER to begin'],'center','center',[0 0 0]);
    
    Screen('Flip',window);
    
    presstocontinue(enterKey,escapeKey);
    
    exploration_satisfied = false;

else
     DrawFormattedText(window,['You scored ' num2str(quiz_accuracy) '% on the quiz',...
        '\n\n Looks like you´ve learnt the maze!',...
       '\n\n Now we´re ready for the next phase of the experiment.',...
       '\n\n\n Press SPACEBAR to continue'],'center','center',[0 0 0]);
   
   Screen('Flip',window);
    
   presstocontinue(spaceKey,escapeKey);
   
   exploration_satisfied = true;
    
end
