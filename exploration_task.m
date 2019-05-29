%% Exploration

DrawFormattedText(window,...
    ['In this experiment, you will play the role of a mouse inside a house.',...
     '\n\n First, you must learn the layout of the rooms in the house. Each room is represented by a unique image.',...
    '\n\n You can use arrow keys to navigate up, left, down, and right.',...
    '\n\n\n\n Press ENTER to begin (or S to skip).'],...
    'center',yCentre-150);
Screen('DrawTexture', window, character_img(1), [], CenterRectOnPointd([0 0 stimulus_size], xCentre, yCentre+200));
Screen('Flip',window);

response = false;
skip = false;
while ~response
    [keyIsDown,secs,keyCode] = KbCheck;
    if keyCode(escapeKey)
        sca
        return
    elseif keyCode(enterKey)
        response = true;
    elseif keyCode(KbName('s'))
        skip = true;
    end
end

if ~skip
    baseRect = CenterRectOnPointd([0 0 stimulus_size], xCentre, yCentre);

    room = start_state;
    end_exploration = false;
    explore_log = [];
    while ~end_exploration

        % Display room & text
        if isempty(explore_log)
            DrawFormattedText(window,'This is your starting room','center',yCentre-200);
        elseif transitioned
            DrawFormattedText(window,'You have entered another room','center',yCentre-200);
        elseif ~transitioned
            DrawFormattedText(window,'Looks like there is a wall here...','center',yCentre-200);
        end

        Screen('DrawTexture', window, img(states(room)), [], baseRect);
        DrawFormattedText(window,['Press the up/down/left/right arrows to move to another room',...
                                 '\n\n\n\n (press ESCAPE to end exploration)'],...
                                 'center',yCentre+200);
        Screen('Flip',window);

        response = false;
        while ~response

            [keyIsDown,secs,keyCode] = KbCheck;

            if keyCode(upKey) || keyCode(downKey) || keyCode(leftKey) || keyCode(rightKey)

                response = true;
                this_action = find(key_map == find(keyCode));

                [transitioned, next_room] = state_transition(this_action,room,transitions,params,1);

            elseif keyCode(escapeKey)
                WaitSecs(1);
                response = true;
                end_exploration = true;
            end
        end

        if ~end_exploration
            explore_log = [explore_log; room, this_action, next_room];
            room = next_room;
        end

        if nExplore ~= Inf
            if size(explore_log,1) == nExplore
                WaitSecs(1);
                end_exploration = true;
            end
        end

    end
end