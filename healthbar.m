% Draw the health status in top left

top_left_char = CenterRectOnPointd([0 0 150 150], 100, 125); % draw the little mouse in the top left
Screen('DrawTexture', window, character_img, [], top_left_char);
fill_length = 150*health;
bar_start = 100 - 150/2;
health_fill = CenterRectOnPointd([0 0 fill_length 30], fill_length/2 + bar_start, 225); % draw the health bar fill
health_rect = CenterRectOnPointd([0 0 150 30], 100, 225); % draw the health bar frame
DrawFormattedText(window,[num2str(round(health*100)) '%'],180,235); % state health in %
DrawFormattedText(window,['Day ' num2str(this_day)],30,40,[0 0 0]); % Day Number
Screen('FillRect', window, [1 0 0], health_fill);
Screen('FrameRect', window, [0 0 0], health_rect, 3);
