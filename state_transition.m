function [transitioned, state2] = state_transition(action,state1,transitions,params)
% Visually present the transitions between two states in the map.

% REQUIRED VARIABLES:

% state1 = starting state (numeric)

% action = 1, 2, 3 or 4 (numeric), indicating the chosen action

% transitions = a list of where every action (column 2 = action) in each state
%               (column 1 = start state) would take you (column 3 = finish state)

% params = a list of extra stuff needed to run this function

img = params.img;
stimulus_size = params.stimsize;
xCentre = params.xCentre;
yCentre = params.yCentre;
window = params.window;
states = params.states;
ifi = params.ifi;
screenXpixels = params.screenX;
screenYpixels = params.screenY;
response_text = params.response_text;
action_labels = params.action_labels;

%%

state2 = transitions(transitions(:,1) == state1 & transitions(:,2) == action,3);
transitioned = state1 ~= state2;

% Draw movement animation
baseRect = CenterRectOnPointd([0 0 stimulus_size], xCentre, yCentre);
Screen('DrawTexture', window, img(states(state1)), [], baseRect); % draw picture again plus chosen direction
DrawFormattedText(window,[response_text{action} ' ' action_labels{action} ' ' response_text{action}],'center',yCentre-200,[0 0 0]); % arrow response
Screen('Flip',window);
WaitSecs(.5);

% with transition = full movement
if any(action == [1 2]) && transitioned % up or down
    
    amplitude = (yCentre + stimulus_size(1))/2;
    swipe_dur = 1; % how many seconds it takes for the stimulus to move completely
    
elseif any(action == [3 4]) && transitioned % left or right
    
    amplitude = (xCentre + stimulus_size(2))/2;
    swipe_dur = 1;
    
    % without transition = wiggle
elseif any(action == [1 2]) && ~transitioned % up or down
    
    amplitude = screenYpixels/16; % just a small movement
    swipe_dur = 1;
    
elseif any(action == [3 4]) && ~transitioned % left or right
    
    amplitude = screenYpixels/16; % just a small movement
    swipe_dur = 1;
    
end

if transitioned
    frequency = 1/swipe_dur/2; % whole duration should do a half cycle
elseif ~transitioned
    frequency = 1/swipe_dur/.25; % whole duration should do 4 full cycles
end
t_idx = 0:ifi:swipe_dur;
angFreq = 2 * pi * frequency;

if ~transitioned
    exp_amp = [];
    for t = 1:length(t_idx)
        exp_amp(t) = amplitude * (.01 ^ -((length(t_idx) - t) / length(t_idx)));
    end
    exp_amp = exp_amp/(max(exp_amp)/amplitude);
end

state1pos = [];
for t = 1:length(t_idx)
    if transitioned
        this_amplitude = amplitude; % amplitude stays constant for transitions (i.e. moves off screen)
        state1pos(t) = this_amplitude * cos(angFreq*t_idx(t)); % starts in trough of cosine wave
    elseif ~transitioned
        this_amplitude = exp_amp(t); % amplitude gets smaller and smaller for 'bounce' effect
        state1pos(t) = this_amplitude * sin(angFreq*t_idx(t)); % starts at zero point of sine wave
    end
end

if transitioned
    if action == 1 % up
        state1pos = fliplr(state1pos);
        state1pos = state1pos + max(state1pos) + yCentre;
        nextstate1pos = state1pos - yCentre - stimulus_size(1);
    elseif action == 2 % down
        state1pos = state1pos - max(state1pos) + yCentre; % 0 to max screen size
        nextstate1pos = state1pos + yCentre + stimulus_size(1); % same, but shifted to the left
    elseif action == 3 % left
        state1pos = fliplr(state1pos);
        state1pos = state1pos + max(state1pos) + xCentre;
        nextstate1pos = state1pos - xCentre - stimulus_size(2);
    elseif action == 4 % right
        state1pos = state1pos + (xCentre-max(state1pos));
        nextstate1pos = state1pos + xCentre + stimulus_size(2);
    end
elseif ~transitioned
    if action == 1 || action == 3 % reverse for down or left
        state1pos = -state1pos;
    end
end

f = 1;
WaitSecs(.5);
while f < length(state1pos)
    
    if transitioned
        
        if any(action == [1 2]) % up/down
            x1 = xCentre;
            x2 = xCentre;
            y1 = state1pos(f);
            y2 = nextstate1pos(f);
        elseif any(action == [3 4]) % left/right
            x1 = state1pos(f);
            x2 = nextstate1pos(f);
            y2 = yCentre;
            y1 = yCentre;
        end
        
        movingRect_state1 = CenterRectOnPointd([0 0 stimulus_size], x1, y1);
        movingRect_nextstate1 = CenterRectOnPointd([0 0 stimulus_size], x2, y2);
        Screen('DrawTexture', window, img(states(state1)), [], movingRect_state1);
        Screen('DrawTexture', window, img(states(state2)), [], movingRect_nextstate1);
        
    elseif ~transitioned
        
        if any(action == [1 2]) % up/down
            x = xCentre;
            y = state1pos(f)+yCentre;
        elseif any(action == [3 4]) % left/right
            x = state1pos(f)+xCentre;
            y = yCentre;
        end
        
        movingRect_state1 = CenterRectOnPointd([0 0 stimulus_size], x, y);
        Screen('DrawTexture', window, img(states(state1)), [], movingRect_state1);
        
    end
    
    vbl  = Screen('Flip', window);
    
    f = f + 1;
    
end

end