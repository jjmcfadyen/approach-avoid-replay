function [training,stateNames,plotdata] = getFL(subject,fl,event)

% directories
dir_meg = fl.dir_processed;

% get neighbours for interpolation
cfg = [];
cfg.method = 'template';
cfg.template = 'ctf275_neighb.mat';
neighbours = ft_prepare_neighbours(cfg);

training = [];
plotdata = cell(1,6); % per condition
for f = 1:length(fl.data)
    
    % load
    load(fullfile(dir_meg,[subject '_FL_' num2str(f) '_osl.mat'])); % loads 'D'
    minfo = load(fullfile(dir_meg,[subject '_FL_' num2str(f) '_osl_info.mat'])); % loads 'dinfo' but change it to 'minfo' 
    minfo = minfo.dinfo;
    
    % interpolate bad sensors
    if ~isempty(minfo.badchannels)
        cfg = [];
        cfg.badchannel = D.label(minfo.badchannels);
        cfg.neighbours = neighbours;
        cfg.senstype = 'meg';
        cfg.grad = D.grad;
        D = ft_channelrepair(cfg,D);
    end
    
    % find MEG sensors (labels that start with 'M')
    sensIdx = startsWith(D.label,'M','IgnoreCase',false);
    
    cfg = [];
    cfg.channel = 'MEG';
    cfg.showcallinfo = 'no';
    P = ft_selectdata(cfg,D);
    
    % get data per trial
    nTrls = length(D.trial);
    nTime = D.time{1}; % must be the same for all trials
    data = nan(nTrls,sum(sensIdx),length(nTime));
    for trl = 1:nTrls
        data(trl,:,:) = D.trial{trl}(sensIdx,:);
    end

    % get label per trial
    L = nan(nTrls,1);
    for trl = 1:nTrls
        idx = zeros(length(minfo.events{trl}),1);
        for i = 1:length(idx)
            if strcmp(minfo.events{trl}(i).type,'image')
                idx(i,1) = 1;
            end
        end
        L(trl,1) = minfo.events{trl}(find(idx)).value;
    end
    
    % check correspondance between MEG and behaviour
    tmp = mean(L == event.table.ID(event.table.Block == f));
    if tmp < 1
        warning('-----------------------------------------------------------------')
        warning(['MEG triggers only match up with ' num2str(round(tmp*100,2)) ' % of behavioural log']);
        warning('-----------------------------------------------------------------')
    end
    
    % remove bad trials
    zrt = abs(zscore(event.table.RT(event.table.Block == f))); % z-scored rt
    minfo.badtrials = [minfo.badtrials find(zrt>3)'];
    if ~isempty(minfo.badtrials)
        data(minfo.badtrials,:,:) = [];
        L(minfo.badtrials,:,:) = [];
        
        cfg = [];
        cfg.trials = setdiff(1:length(D.trial),minfo.badtrials);
        cfg.showcallinfo = 'no';
        P = ft_selectdata(cfg,P);
    end
    
    % save to variable
    if f == 1
        training.x = nTime;
        training.D = data;
        training.L = L;
        for st = 1:6
            cfg = [];
            cfg.trials = find(L == st);
            cfg.showcallinfo = 'no';
            plotdata{st} = ft_selectdata(cfg,P);
        end
    else
        training.D = [training.D; data];
        training.L = [training.L; L];
        for st = 1:6
            cfg = [];
            cfg.trials = find(L == st);
            cfg.showcallinfo = 'no';
            plotdata{st} = ft_appenddata([],plotdata{st},ft_selectdata(cfg,P));
        end
    end
end

stateNames = {};
for st = 1:6
    stateNames{st} = event.table.Image{find(event.table.ID == st,1,'first')};
end

% % Plot
%{

avgs = cell(1,6);
for st = 1:6
    avgs{st} = ft_timelockanalysis([],plotdata{st});
end

cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.showlabels = 'yes';
cfg.linecolor = colours(6,'rainbow');
figure; ft_multiplotER(cfg,avgs{:});

%}

end