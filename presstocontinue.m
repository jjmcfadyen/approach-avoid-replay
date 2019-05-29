function presstocontinue(continue_key,abort_key)

% continue_key = the code for the key you're interested in people pressing to continue
%                Should be whatever KbName('s') gives you
% abort_key = same as above, but used for exiting the experiment entirely

response = false;
    
[~,~,keyCode] = KbCheck;
while ~response
    old_keyCode = keyCode;
    [~,~,keyCode] = KbCheck;
    keyDiff = keyCode - old_keyCode; % 1 = newly pressed, -1 = newly released
    if keyDiff(abort_key)
        sca
        return
    elseif keyDiff(continue_key) == -1 % looks for a RELEASE
        response = true;
    end
end

end