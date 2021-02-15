function [data,badidx] = getTask(subject,fl,event,epoch)

% directories
dir_meg = fl.dir_processed;

% parameters
responses = [4 1]; % 4 = risky, 1 = safe

% get neighbours for interpolation
cfg = [];
cfg.method = 'template';
cfg.template = 'ctf275_neighb.mat';
cfg.feedback = 'no';
warning off
neighbours = ft_prepare_neighbours(cfg);
warning on

data = [];
badidx = [];
for f = 1:length(fl.data)
    
    % load
    load(fullfile(dir_meg,[subject '_Task_' num2str(f) '_' epoch '.mat'])); % loads 'D'
    warning off
    minfo = load(fullfile(dir_meg,[subject '_Task_' num2str(f) '_' epoch '_info.mat'])); % loads 'dinfo' but change it to 'minfo' 
    warning on
    minfo = minfo.dinfo;
    
    % interpolate bad sensors
    if ~isempty(minfo.badchannels) && length(unique(D.grad.chanpos)) > 1 % if the channel positions are all zeros, can't do interpolation
        cfg = [];
        cfg.badchannel = D.label(minfo.badchannels);
        cfg.neighbours = neighbours;
        cfg.senstype = 'meg';
        cfg.grad = D.grad;
        warning off
        D = ft_channelrepair(cfg,D);
        warning on
    end
    
    % find MEG sensors (labels that start with 'M')
    sensIdx = startsWith(D.label,'M','IgnoreCase',false);
    
    cfg = [];
    cfg.channel = 'MEG';
    cfg.showcallinfo = 'no';
    P = ft_selectdata(cfg,D);
    
    % check for doubled-up samples
    [u,idx,uidx] = unique(D.sampleinfo,'rows');
    D.time = D.time(idx);
    D.trial = D.trial(idx);
    D.sampleinfo = D.sampleinfo(idx,:);
    D.trialinfo = D.trialinfo(idx,:);
    minfo.badtrials = unique(uidx(minfo.badtrials));
    minfo.events = minfo.events(idx);
    
    % get data & time info per trial
    nTrls = length(D.trial);
    d = cell(1,nTrls);
    timeinfo = cell(1,nTrls);
    for trl = 1:nTrls
        d{trl} = D.trial{trl}(sensIdx,:);
        if size(d{trl},1) ~= nTrls
            d{trl} = d{trl}';
        end
        timeinfo{trl} = D.time{trl};
    end
    
    % get label per trial & crosscheck with 'event' variable & log bad trials
    L = nan(nTrls,2); % first col = exp trial, second col = value
    for trl = 1:nTrls
        ev = minfo.events{trl};
        switch epoch
            case 'decision'
                cev = [];
                for i = 1:length(ev)
                    if contains(ev(i).type,'choice')
                        cev = ev(i);
                        break;
                    end
                end
                rev = [];
                for i = 1:length(ev)
                    if contains(ev(i).type,'response')
                        rev = ev(i);
                        break;
                    end
                end
                L(trl,1) = cev.value; % event.table.ExpTrial number
                if ~isempty(rev)
                    L(trl,2) = rev.value; % response value (see 'responses' variable for mappings)
                else
                    L(trl,2) = NaN;
                end
            case 'transition' % label is event.table.ExpTrial number
                tev = [];
                for i = 1:length(ev)
                    if contains(ev(i).type,'transition')
                        tev = ev(i);
                        break;
                    end
                end
                iev = [];
                for i = 1:length(ev)
                    if contains(ev(i).type,'image')
                        iev = ev(i);
                        break;
                    end
                end
                L(trl,1) = tev.value; % event.table.ExpTrial number
                if ~isempty(iev)
                    L(trl,2) = str2num(iev.type(8)); % path from subsequent image type
                else
                    L(trl,2) = NaN;
                end
            case 'image' % label is the state number (1-6)
                for i = 1:length(ev)
                    if contains(ev(i).type,'image')
                        ev = ev(i);
                        break;
                    end
                end
                pnum = event.table.Transition(event.table.ExpTrial == ev.value);
                snum = str2num(ev.type(end));
                L(trl,1) = ev.value;
                if pnum == 1
                    L(trl,2) = snum;
                elseif pnum == 2
                    L(trl,2) = snum+3;
                end
            case 'outcome' % label is event.table.ExpTrial number   
                for i = 1:length(ev)
                    if contains(ev(i).type,'outcome')
                        ev = ev(i);
                        break;
                    end
                end
                L(trl,1) = ev.value; % label is event.table.ExpTrial number
        end
    end
    
    % manual error fix
    if strcmp(epoch,'decision')
        if strcmp(subject,'680913')
            idx = ismember(L(:,1),[68 122 140]);
            L(idx,2) = responses(1); % manual error fix
        elseif strcmp(subject,'396430')
            idx = ismember(L(:,1),[3, 50]);
            L(idx,2) = responses(1); % manual error fix (first 'safe' response not registered) 
        end
    end
    
    % check against behavioural file
    v = true;
    this_T = event.table(ismember(event.table.ExpTrial,L(:,1)),:);
    switch epoch
        case 'decision'
            tmp = L(:,2);
            for r = 1:length(responses)
                L(tmp==responses(r),2) = r;
            end
            if any((this_T.Choice(~isnan(L(:,2))) == L(~isnan(L(:,2)),2)) == 0)
                v = false;
            end
        case 'transition'
            if any(this_T.Choice == 2) % transitions only happen for risky (1) not safe (2) choices
                v = false;
            end
            if any((this_T.Transition(~isnan(L(:,2))) == L(~isnan(L(:,2)),2)) == 0)
                v = false;
            end
        case 'image'
            uI = unique(L(:,1));
            for i = 1:length(uI)
                tmp = this_T.Transition(this_T.ExpTrial == uI(i));
                if (tmp==1 && ~any(ismember(L(L(:,1) == uI(i),2),1:3))) || ...
                   (tmp==2 && ~any(ismember(L(L(:,1) == uI(i),2),4:6)))  
                    v = false;
                end
            end
    end
    if ~v
        error('Something went wrong when aligning the MEG trials with behaviour.')
    end
    
    % log bad trials
    bt = zeros(nTrls,1);
    bt(minfo.badtrials,1) = 1;
    
    % save to variable
    if f == 1
        data.x = timeinfo;
        data.D = d;
        data.L = L;
        try
            data.Fs = D.fsample;
        catch
            warning(['Could not get sampling rate for ' subject ', file ' num2str(f)])
        end
        badidx = bt;
    else
        data.x = [data.x, timeinfo];
        data.D = [data.D, d];
        data.L = [data.L; L];
        badidx = [badidx; bt];
        if ~isfield(data,'Fs')
            try
                data.Fs = D.fsample;
            catch
                warning(['Could not get sampling rate for ' subject ', file ' num2str(f)])
            end
        end
    end
end

end