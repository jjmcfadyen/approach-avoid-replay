% Main preprocessing script

clear all
clc

%% Directories

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_meg = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';

cd D:\2020_RiskyReplay\approach-avoid-replay

%% Parameters

addpath('utils');
addpath('preprocessing')

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

% Downsampling
Fs = [100 600];

%% Convert data to SPM

for s = 1:N
    
    % Get list of raw data files
    cd(fullfile(dir_raw,subjects{s}))
    idx = find(parameters.subjectID==str2double(subjects{s}));
    filelist = parameters.rawfile(idx);
    
    % Set output directory
    thisoutput = fullfile(dir_meg,'1_converted',subjects{s});
    if ~exist(thisoutput)
        mkdir(thisoutput)
    end
    
    % Convert each file
    for f = 1:length(filelist)
        
        disp('===========================')
        disp(['CONVERTING ' subjects{s} ', run ' num2str(f) ' of ' num2str(length(filelist))])
        disp('===========================')
        S = [];
        S.dataset = fullfile(dir_raw,subjects{s},filelist{f});
        D = spm_eeg_convert(S);
        D.move(fullfile(thisoutput,...
            ['spm_' subjects{s} '_' parameters.task{idx(f)} '_r' num2str(parameters.block(idx(f))) '.mat']))
    end
end

%% Identify triggers in photodiode and crop data

for s = 1:N
    
    % Get list of raw data files
    idx = find(parameters.subjectID==str2double(subjects{s}));
    filelist = parameters.rawfile(idx);
    
    % Get behavioural logs
    behav = [];
    behav.FL = readtable(fullfile(dir_behav,subjects{s},[num2str(str2double(subjects{s})) '_fl.csv']));
    behav.FL.Subject = repmat(subjects(s),size(behav.FL,1),1);
    behav.task = parse_behav(subjects{s},dir_behav);
    behav.task.Subject = cellstr(behav.task.Subject);
    
    save(fullfile(dir_behav,subjects{s},[subjects{s} '_parsedBehav.mat']),'behav','parameters');
    
    % Set input directory
    thisinput = fullfile(dir_meg,'1_converted',subjects{s});
    if ~exist(thisinput)
        mkdir(thisinput)
    end
    
    % Set output directory
    thisoutput = fullfile(dir_meg,'2_cropped',subjects{s});
    if ~exist(thisoutput)
        mkdir(thisoutput)
    end
    
    % Get triggers from each file
    for f = 1:length(filelist)
        
        disp('==========================================')
        disp(['GETTING TRIGGERS FOR ' subjects{s} ', run ' num2str(f) ' of ' num2str(length(filelist))])
        disp('==========================================')
        
        thistask = parameters.task{idx(f)};
        thisblock = parameters.block(idx(f));
        
        % Load converted data
        D = spm_eeg_load(fullfile(thisinput,...
            ['spm_' subjects{s} '_' thistask '_r' num2str(thisblock) '.mat']));
        
        % Select relevant behavioural log
        thisbehav = [];
        switch thistask
            case 'FL'
                thisbehav = behav.FL(behav.FL.Block==thisblock,:);
            case 'task'
                if thisblock==0
                    thisbehav = behav.task(behav.task.Practice==1,:);
                else
                    thisbehav = behav.task(behav.task.Practice==0 & behav.task.Block==thisblock,:);
                end
        end
        
        % Get triggers
        [triggers,croptime] = parse_triggers(D,thisbehav,thistask);
        
        % Crop the MEG data
        S = [];
        S.D = D;
        S.timewin = [0 croptime*1000]; % in ms
        cropped = spm_eeg_crop(S);
        
        % Insert the triggers as events and save
        cropped = insert_events(cropped,triggers,thistask);
        cropped.move(fullfile(thisoutput,...
            ['cropped_' subjects{s} '_' thistask '_r' num2str(thisblock) '.mat']));
        
        save(fullfile(thisoutput,...
            ['triggers_' subjects{s} '_' thistask '_r' num2str(thisblock) '.mat']),'triggers');
        
    end
end

%% Visually validate the EOG channels

loadEOG = false; % if EOG channels have already been saved, load them instead

if ~loadEOG
    figure
    set(gcf,'position',[5 219 1431 777])
end
for s = [9 20]
   
    dir_input = fullfile(dir_meg,'2_cropped',subjects{s});
    filelist = dir(fullfile(dir_input,'cropped*.mat'));
    
    for f = 1:length(filelist)
        
        [subject, task, run] = split_filename(filelist(f).name);
    
        D = spm_eeg_load(fullfile(filelist(f).folder,filelist(f).name));

        disp([subject ' ' task ' ' num2str(run) ' sum = ' num2str(sum(contains(D.chantype,'EOG')))])
        
        if f==1 && ~loadEOG
            clf
            uadchans = find(contains(D.chanlabels,'UADC'));
            pdata = squeeze(D(uadchans,:,:));
            nChan = size(pdata,1);
            for i = 1:nChan
                if nChan+1 <= 4
                   subplot(size(pdata,1)+1,1,i)
                else
                    subplot(round((size(pdata,1)+1)/2),2,i)
                end
                plot(D.time,pdata(i,:))
                title(D.chanlabels(uadchans(i)))
            end
            imagesc(corr(pdata')); colormap('hot'); caxis([0 1])
            drawnow

            eog = inputdlg('EOG Channels:','EOG',1,{'UADC001 UADC002 UADC003'});
            if length(eog)==1
                eog = strsplit(eog{1},' ');
            end

            save(fullfile(dir_input,[subjects{s} '_eogLabels.mat']),'eog');
        elseif loadEOG
            load(fullfile(dir_input,[subjects{s} '_eogLabels.mat']));
        end

        D = chantype(D,find(contains(D.chanlabels,eog)),'EOG');
        D.save;
    end
end

%% Generate jobs for OSL on cluster (filter/downsample/ICA)

addpath('D:\2020_RiskyReplay\approach-avoid-replay\preprocessing\batch')

for s = [9 20]
   
    filelist = dir(fullfile(dir_meg,'2_cropped',subjects{s},'cropped*.mat'));
%     fileidx = find(parameters.subjectID==str2double(subjects{s}));
%     filelist = struct();
%     for f = 1:length(fileidx)
%         filelist(f).name = ['spm_' subjects{s} '_' parameters.task{f} '_r' num2str(parameters.block(fileidx(f))) '.mat'];
%     end

    for f = 1:length(filelist)
        generate_jobs_preprocess(filelist(f).name);
    end
    
end

%% Epoch the data

planningType = 'during'; % 'during' planning (main analysis) or 'post' planning (transition to outcome)

% for plotting merged conditions in Fieldtrip
addpath('D:\Toolboxes\fieldtrip-20191119\template\layout')
addpath('D:\Toolboxes\fieldtrip-20191119\template\neighbours')

cfg = [];
cfg.method = 'template';
cfg.layout = 'CTF275.lay';
neighbours = ft_prepare_neighbours(cfg);

trigger_errors = [];
for ds = 1%:length(Fs)
    for s = [1:2 4:N]%1:N % still need to do 20 for 600 Hz

        disp('==========================================================================================')
        disp('==========================================================================================')
        disp('==========================================================================================')
        disp(['=========== ' subjects{s} '========================================================================'])
        disp('==========================================================================================')
        disp('==========================================================================================')
        disp('==========================================================================================')
        
        %% Epoch
        dir_output = fullfile(dir_meg,['6_epoched_ds-' num2str(Fs(ds)) 'Hz'],subjects{s});
        if ~exist(dir_output)
            mkdir(dir_output)
        end
        
        % Get behavioural logs for main task
        behav = parse_behav(subjects{s},dir_behav);
        
        % Get files
        filelist = dir(fullfile(dir_meg,['5_ICA_ds-' num2str(Fs(ds)) 'Hz'],subjects{s},'*.mat'));
        if isempty(filelist)
            error(['No files found for ' subjects{s} ' (' num2str(Fs(ds)) ' Hz)'])
        end
        
        for f = 5:length(filelist)
            
            [subject, task, run] = split_filename(filelist(f).name);

            % Get events from data file
            D = spm_eeg_load(fullfile(filelist(f).folder,filelist(f).name));
            events = D.events;

            vals = extractfield(events,'value');
            if ~isnumeric(vals)
                tmp = nan(length(vals),1);
                tmp(cellfun(@isnumeric,vals)) = cell2mat(vals(cellfun(@isnumeric,vals)));
                vals = tmp;
            end
            if size(vals,2)>size(vals,1)
                vals = vals';
            end
            
            types = extractfield(events,'type')';
            onsets = extractfield(events,'time')';
            
            uVals = unique(vals(~isnan(vals)));

            % Define trials
            trl = [];
            switch task
                case 'FL'
                    timewin = [-.1 .5];
                    for st = 1:length(uVals)
                        idx = strcmp(types,'picture') & vals==uVals(st);
                        trl = [trl; D.indsample(onsets(idx)+timewin(1))',... % onset minus baseline
                                    D.indsample(onsets(idx)+timewin(2))',... % onset plus epoch length
                                    repmat(timewin(1)*D.fsample,sum(idx),1),... % offset of zero point
                                    repmat(uVals(st),sum(idx),1)]; % condition number
                    end
                case 'task'
                    for st = 1:length(uVals)
                        
                        thisval = sprintf('%04d',uVals(st));
                        thisblock = str2double(thisval(1:2));
                        thistrial = str2double(thisval(3:end));
                        
                        if thisblock==0
                            thisidx = behav.Practice==1 & behav.Block==1 & behav.Trial==thistrial;
                        else
                            thisidx = behav.Practice==0 & behav.Block==thisblock & behav.Trial==thistrial;
                        end
                        switch planningType
                            case 'during'
                                
                                editedtrigger = false;
                                idx = strcmp(types,'decision') & vals==uVals(st);
                                
                                thisrt = behav.RT(thisidx);
                                if find(idx) < length(vals)
                                    ii = find(idx)+1;
                                    while strcmp(types{ii},'artefact_OSL')
                                        ii = ii+1;
                                        if ii > length(vals)
                                            ii = NaN;
                                            break
                                        end
                                        if strcmp(types{ii},'decision')
                                            ii = NaN;
                                            break
                                        end
                                    end
                                    if ~isnan(ii)
                                        if events(ii).time - events(find(idx)).time < thisrt
                                            thisrt = events(find(idx)+1).time - events(find(idx)).time - 1;
                                            editedtrigger = true;
                                            warning('RT is longer than the time between trial start and transition/outcome')
                                        end
                                    end
                                end
                                
                                if Fs(ds) == 100
                                    bc = [-.5 .5]; % 500 ms before trial onset to 500 ms after response
                                elseif Fs(ds) == 600
                                    bc = [-1 1]; % 1000 ms before trial onset to 1000 ms after response
                                end
                                timewin = [bc(1) thisrt+bc(2)];
                                
                                if abs((timewin(end)-bc(2)) - thisrt) > .5
                                    error('RT does not match up')
                                end
                                
                                trigger_errors = [trigger_errors; ...
                                    array2table([s f st timewin(end)-bc(2) behav.RT(thisidx) editedtrigger],'variablenames',{'Subject','File','Event','TriggerWindow','RT','Edited'})];
                                
                                thistrl = [D.indsample(onsets(idx)+timewin(1))',... % onset minus baseline
                                    D.indsample(onsets(idx)+timewin(2))',... % onset plus epoch length
                                    repmat(timewin(1)*D.fsample,sum(idx),1),... % offset of zero point
                                    repmat(thisblock,sum(idx),1),... % block
                                    repmat(thistrial,sum(idx),1)]; % trial;
                                
                                negidx = isnan(thistrl(1));
                                if ~negidx
                                    negidx = D.time(thistrl(1)) < 0;
                                end
                                if negidx
                                    timewin(negidx) = 1;
                                end
                                
                                posidx = isnan(thistrl(2));
                                if ~posidx
                                    posidx = thistrl(2) > length(D.time);
                                end
                                if posidx
                                    timewin(posidx) = length(D.time);
                                end
                                
                                thistrl = [D.indsample(onsets(idx)+timewin(1))',... % onset minus baseline
                                    D.indsample(onsets(idx)+timewin(2))',... % onset plus epoch length
                                    repmat(timewin(1)*D.fsample,sum(idx),1),... % offset of zero point
                                    repmat(thisblock,sum(idx),1),... % block
                                    repmat(thistrial,sum(idx),1)]; % trial;
                                
                                trl = [trl; thistrl];
                                
                            case 'post'
                                
                                % get the first & last trigger AFTER the decision for this trial
                                allthisval = find(vals==uVals(st)); % every trigger for this trial
                                
                                if length(allthisval) > 1
                                    thisonset = extractfield(events,'time');
                                    thisonset = thisonset([allthisval(2) allthisval(end)]); % take first trigger after decision is made to last one for this trial
                                    thisonset(end) = thisonset(end)+2; % add 2 seconds after outcome presentation

                                    % create an index so that it's timelocked to the FIRST trigger after decision
                                    idx = strcmp(types,types{allthisval(2)}) & vals==uVals(st);
                                    if sum(idx)>1
                                        error('Error with idx - too many post-decision triggers to choose from')
                                    end
                                    timewin = [0 thisonset(end)-thisonset(1)];

                                    trl = [trl; D.indsample(thisonset(1))',... % onset
                                            D.indsample(thisonset(end))',... % offset
                                            0,... % offset of zero point
                                            thisblock,... % block
                                            thistrial]; % trial
                                end
                        end
                    end
            end
            
            if ~isempty(trl)
                
                trl = sortrows(trl,1);
                condlist = trl(:,4:end);
                trl = trl(:,1:3);

                % Convert to fieldtrip
                ft = ftraw(D);

                % Epoch using 'trl' variable
                cfg = [];
                cfg.trl = trl;

                epoched = ft_redefinetrial(cfg,ft);
                epoched.trialinfo = condlist;

                % Only select MEG channels
                cfg = [];
                cfg.channel = D.chanlabels(find(contains(D.chantype,'MEGGRAD')));

                epoched = ft_selectdata(cfg,epoched);

                if isfield(epoched,'elec')
                    epoched = rmfield(epoched,'elec'); % subject 707132 had random channels included that were identified as EEG channels
                end

                % Interpolate bad channels
                cfg = [];
                cfg.badchannel = D.chanlabels(D.badchannels)';
                cfg.neighbours = neighbours;

                epoched = ft_channelrepair(cfg,epoched);

                % Save
                plantitle = '';
                if strcmp(planningType,'post')
                    plantitle = 'post-';
                end
                save(fullfile(dir_output,...
                    ['epoched_' num2str(Fs(ds)) 'Hz_' subjects{s} '_' plantitle task '_r' num2str(run) '.mat']),'epoched');
            else
                warning(['No events found for ' subjects{s} ', ' filelist(f).name])
            end
        end
        
        clear D
        clear epoched
        clear ft
        
        %% Merge runs
        
        dir_output = fullfile(dir_meg,['7_merged_ds-' num2str(Fs(ds)) 'Hz']);
        if ~exist(dir_output)
            mkdir(dir_output)
        end
        
        tasks = {'FL','task'};
        
        for t = 2%1:length(tasks)
            
            plantitle = '';
            if strcmp(planningType,'post') && strcmp(tasks{t},'task')
                plantitle = 'post-';
            end
            
            filelist = dir(fullfile(dir_meg,['6_epoched_ds-' num2str(Fs(ds)) 'Hz'],subjects{s},...
                ['epoched_' num2str(Fs(ds)) 'Hz_' subjects{s} '_' plantitle tasks{t} '*.mat']));
    
            % merge
            tomerge = cell(1,length(filelist));
            for f = 1:length(filelist)
                tmp = load(fullfile(filelist(f).folder,filelist(f).name));
                tomerge{f} = tmp.epoched;
            end
            
            merged = ft_appenddata([],tomerge{:});
            
            % save file
            save(fullfile(dir_output,[subjects{s} '_' plantitle tasks{t} '_' num2str(Fs(ds)) 'Hz.mat']),'merged','-v7.3');
            
%             % if main task, make a response-locked version as well
%             if t==2 && strcmp(planningType,'during')
%                
%                 cfg = [];
%                 cfg.offset = nan(length(merged.trial),1);
%                 for trl = 1:length(cfg.offset)
%                     cfg.offset(trl,1) = -findMin(merged.time{trl}(end)-1,merged.time{trl}); 
%                 end
%                 
%                 resplocked = ft_redefinetrial(cfg,merged);
%                 
%                 merged = resplocked;
%                 save(fullfile(dir_output,[subjects{s} '_response_' num2str(Fs(ds)) 'Hz.mat']),'merged','-v7.3');
%                 
%             end
            
        end
        
        % plot
%         for t = 2%1:length(tasks)    
% 
%             plantitle = '';
%             if strcmp(planningType,'post') && strcmp(tasks{t},'task')
%                 plantitle = 'post-';
%             end
%             
%             load(fullfile(dir_output,[subjects{s} '_' plantitle tasks{t} '_' num2str(Fs(ds)) 'Hz.mat'])); % loads 'merged' variable
%             
%             if t==1
%                 uCon = unique(merged.trialinfo);
%                 pdata = cell(1,length(uCon));
%                 for c = 1:length(uCon)
%                     cfg = [];
%                     cfg.trials = find(merged.trialinfo==uCon(c));
%                     pdata{c} = ft_timelockanalysis(cfg,merged);
%                 end
%             else
%                 pdata = {ft_timelockanalysis([],merged)};
%             end
% 
%             figure
%             cfg = [];
%             cfg.layout = 'CTF275.lay';
%             cfg.parameter = 'avg';
%             if t==1
%                 cfg.linecolor = colours(length(uCon),'viridis');
%             end
%             ft_multiplotER(cfg,pdata{:})
%             sgtitle(tasks{t})
%             
%         end
    end
end

