%% Task instructions

left_pos = (screenXpixels * .4)-stimulus_size(1);
right_pos = (screenXpixels * .6)+stimulus_size(1);

left_rect = CenterRectOnPointd([0 0 300 300], left_pos, yCentre+150);
right_rect = CenterRectOnPointd([0 0 300 300], right_pos, yCentre+150);
mid_rect = CenterRectOnPointd([0 0 300 300], xCentre, yCentre+150);

understand = false;
    while ~understand

    % Screen 1 - now there are outcomes
    DrawFormattedText(window,['We are now about to start the test phase.',...
                              '\n\n In this phase, you will navigate through the rooms of the house you just learnt.',...
                              '\n\n But this time, we have introduced a reward (CHEESE) and a threat (PREDATOR).',...
                              '\n\n\n Press SPACEBAR to continue'],...
                      'center','center',[0 0 0]);
    Screen('Flip',window);
    presstocontinue(spaceKey,escapeKey);

    % Screen 2 - show the cheese & predator
    DrawFormattedText(window,['You will have to use trial and error to figure out how to get to the CHEESE and avoid the PREDATOR.'],'center',yCentre-150,[0 0 0]);
    Screen('DrawTexture', window, outcome_stimuli(1), [], left_rect);
    Screen('DrawTexture', window, outcome_stimuli(2), [], right_rect);
    DrawFormattedText(window,'Press SPACEBAR to continue.','center',yCentre+150,[0 0 0]);
    Screen('Flip',window);
    presstocontinue(spaceKey,escapeKey);

    % Screen 3 - health bar introduction
    health = .5;
    this_day = 1;
    healthbar; % runs external script
    if loss_mag == Inf
        predator_text = 'die';
    else predator_text 'lose health';
    end
    DrawFormattedText(window,['For 5 days, you must survive in the house.',...
                              '\n\n You will start with a certain amount of health (displayed in the top left)',...
                              '\n\n You can regain health by finding the cheese or getting lucky and finding crumbs.',...
                              '\n\n You will ' predator_text ' if you are attacked by the predator.',...
                              '\n\n If you don’t find any cheese that day, you will go hungry and lose health.',...
                              '\n\n Your goal is to have as much health as possible by the end of the 5 days.',...
                              '\n\n Everything resets each week (i.e. every 5 days).',...
                              '\n\n\n Press SPACEBAR to continue.'],'center','center',[0 0 0]);
    Screen('Flip',window);    
    presstocontinue(spaceKey,escapeKey);

    % Screen 4 - probability introduction
    DrawFormattedText(window,['You must use trial and error to figure out which path leads to the PREDATOR, and which path leads to the CHEESE.',...
                              '\n\n\n On different days, the predator could either be ASLEEP, HALF-AWAKE, or ALERT.',...
                              '\n\n If you walk directly into the predator, it will definitely attack you no matter how awake it is.',...
                              '\n\n If the predator is HALF-AWAKE or ALERT, it is more likely to detect you and attack you from adjacent rooms.',...
                              '\n\n\n Press SPACEBAR to continue.'],'center',yCentre-200,[0 0 0]);
    Screen('DrawTexture', window, threatprob_stimuli(1), [], left_rect);
    Screen('DrawTexture', window, threatprob_stimuli(2), [], mid_rect);
    Screen('DrawTexture', window, threatprob_stimuli(3), [], right_rect);
    Screen('Flip',window);    
    presstocontinue(spaceKey,escapeKey);

    % About to start practice
    DrawFormattedText(window,['At the start of each day, you will be given ' num2str(resp_time) ' seconds to enter in 3 moves.',...
                              '\n\n ** Sometimes, the computer will instruct you to do a specific journey **',...
                              '\n\n Then, you will see your journey and whether you won or lost money.',...
                              '\n\n\n Press ENTER to begin the practice or R to repeat instructions.'],'center','center',[0 0 0]);
    Screen('Flip',window);

    response = false;
    [~,~,keyCode] = KbCheck;
    while ~response
        old_keyCode = keyCode;
        [~,~,keyCode] = KbCheck;
        keyDiff = keyCode - old_keyCode; % 1 = newly pressed, -1 = newly released
        if keyDiff(escapeKey)
            sca
            return
        elseif keyDiff(enterKey) == -1 % looks for a RELEASE
            response = true;
            understand = true;
        elseif keyDiff(KbName('r')) == -1
            response = true;
            understand = false;
        end
    end
    
end
