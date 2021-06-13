function D = insert_events(D,triggers,task)

nTrigs = size(triggers,1);

% Create events structure
ev = struct('type',[],'value',[],'time',[],'duration',[],'offset',[]);
for i = 1:nTrigs
    ev(i).type = triggers.type{i};
    switch task
        case 'FL'
            ev(i).value = triggers.imgNum(i);
        case 'task'
            if all(triggers.practice==1)
                ev(i).value = triggers.trial(i); % practice block, value = trial number
            else
                ev(i).value = triggers.trial(i) + 100*triggers.block(i); % block number * 100 + trial number (e.g., block 2, trial 6 = 206)
            end
    end
    ev(i).time = D.time(triggers.sOnset(i));
    ev(i).duration = triggers.tDur(i);
    ev(i).offset = 0;
end

% Add this to the data
D = events(D,1,ev);
D.save;

end