function preprocessMEG_2(subjects, session)

% subjects = cell array of strings - e.g., {'088674','707132'}
% session = 'FL' or 'Task'

% Uses: output from preprocessMEG_1.m script:
% ---------- [SUBJECT]_[SESSION]_endsamples-[FS]Hz.mat
% ---------- [SUBJECT]_[SESSION]_triggers-[FS]Hz.mat
% Produces: cleaned data (cropped, added events, filtered, downsampled)
% ---------- 
% ---------- 

%% Set up

switch nargin
    case 0
        error('No subjects or session specified')
    case 1
        error('No session specified (FL or Task)');
end

if ~iscell(subjects)
    subjects = {subjects};
end

% directories
cd ..
addpath('utils')
addpath('preprocessing')

[sinfo,dinfo] = dir_cfg;

addpath(genpath(dinfo.tb_osl));
osl_startup();

% Parameters
Fs = 600; % frequency to downsample to (in Hz) - 100 or 600

%% Run

for subj = 1:length(subjects)
   
    subject = subjects{subj};
    [fl,event,in] = pp_cfg(subject,session);
    disp(['----- SUBJECT:  ' subject])
    disp(['----- SESSION:  ' session])
    
    %% Load end samples & triggers
    % (produced by preprocessingMEG_1.m)
    
    fname = dir( fullfile(dinfo.data_meg_pp1,subject,[subject '_' session '_endsamples-*Hz.mat']) );
    if length(fname) > 1
        error(['More than one sampling rate detected for ' subject]);
    end
    load( fullfile(fname.folder,fname.name) );
    
    fname = dir( fullfile(dinfo.data_meg_pp1,subject,[subject '_' session '_triggers-*Hz.mat']) );
    if length(fname) > 1
        error(['More than one sampling rate detected for ' subject]);
    end
    T = load( fullfile(fname.folder,fname.name) );
    T = T.triggerTable;
    
    % Fix label cells
    try
        for i = 1:size(T,1)
            if size(T.label{i},2) == 1
                T.label{i} = T.label{i}{1};
            end
        end
    catch
    end

    % Check for duplicates
    try
        T = unique(T,'rows');
    catch
        idx = unique(T(:,1:3),'rows');
        if size(T,1) ~= size(idx,1)
            error('Duplicate rows in T');
        end
    end
    
    %% Set up directories
    
    dir_output = fullfile(dinfo.data_meg_pp2,subject);
    if ~exist(dir_output)
        mkdir(dir_output);
    end
    
    %% Run each MEG block
    
    tic
    for f = 1:length(fl.data)

        %% Import
        
        outfile = fullfile(dir_output,[subject '_' session '_' num2str(f)]);
        raw = osl_import(fullfile(fl.dir_meg,fl.data{f}),...
                         'outfile',outfile,...
                         'other_channels',in.eog);

        % identify the EOG channels (EyeLink data)
        for i = 1:length(in.eog)
            raw = raw.chantype(find(strcmp(raw.chanlabels,in.eog{i})),'EOG');
            raw = raw.chanlabels(find(strcmp(raw.chanlabels,in.eog{i})),['EOG' num2str(i)]);
        end

        %% Crop

        % crop
        S = [];
        S.D = raw;
        S.timewin = [0 endSamples(f)/raw.fsample*1000]; % in ms
        D = spm_eeg_crop(S); % appends 'p'

        %% Create events
        
        % Crop for just this file
        if strcmp(session,'FL')
            idx = find(event.table.Block == fl.block_idx(f));
        else
            idx = find(event.table.Practice == fl.block_idx(f,2) & event.table.Block == fl.block_idx(f,1));
        end
        triggerTable = T(T.trial >= idx(1) & ...
                         T.trial <= idx(end),:);

        if isempty(triggerTable)
            error('empty trigger table')
        end
                     
        % decide on the events
        if strcmp(session,'FL')

            % get state names
            stateNames = cell(length(unique(event.table.Image)),1);
            for st = 1:length(stateNames)
                stateNames{st,1} = event.table.Image{find(event.table.ID == st,1,'first')};
            end

            % get index of actual image presentations
            idx = zeros(size(triggerTable,1),1);
            for st = 1:length(stateNames)
                idx(strcmp(triggerTable.label,stateNames{st}),1) = 1;
            end
            triggerTable = triggerTable(find(idx),:);

            ev = struct('type',[],'value',[],'time',[],'duration',[],'offset',[]);
            for trl = 1:size(triggerTable,1)
                ev(trl).type = 'image';
                ev(trl).value = find(strcmp(stateNames,triggerTable.label{trl})); % 1 (risky) or 2 (safe)
                ev(trl).time = D.time(triggerTable.onset(trl));
                ev(trl).duration = 0.05;
                ev(trl).offset = 0;
            end

        elseif strcmp(session,'Task')

            % -- GET CHOICES --
            
            choiceIdx = strcmp(triggerTable.label,'Risky choice') | strcmp(triggerTable.label,'Safe choice');
            transIdx = strcmp(triggerTable.label,'Door 1') | strcmp(triggerTable.label,'Door 2') | ...
                              strcmp(triggerTable.label,'Supply room');
                          
            if sum(choiceIdx) == 0
                error('no choices made')
            end

            % look for choices where there's an outcome trigger straight after
            choiceIdx = [0; choiceIdx(1:end-1)] == transIdx;
            choiceTable = triggerTable(choiceIdx,:);

            choiceNames = {'Risky choice','Safe choice'};

            ev = struct('type',[],'value',[],'time',[],'duration',[],'offset',[]);
            cc = 0;
            for trl = 1:size(choiceTable,1)
                if any(strcmp(choiceNames,choiceTable.label(trl)))
                    cc = cc + 1;
                    ev(cc).type = 'choice';
                    ev(cc).value = choiceTable.trial(trl); % trial number in experiment
                    ev(cc).time = D.time(choiceTable.onset(trl));
                    ev(cc).duration = 0.05;
                    ev(cc).offset = 0;
                end
            end

            % -- GET TRANSITION --
            transIdx = [0; choiceIdx(1:end-1)] & ....
                strcmp(triggerTable.label,'Door 1') | strcmp(triggerTable.label,'Door 2') | strcmp(triggerTable.label,'Supply room'); % only get transitions after valid choices
            transTable = triggerTable(transIdx,:);
            
            transNames = {'Door 1','Door 2','Supply Room'};
            
            for trl = 1:size(transTable,1)
                if any(strcmp(transNames,transTable.label(trl)))
                    cc = cc + 1;
                    ev(cc).type = 'transition';
                    ev(cc).value = transTable.trial(trl); % trial number in experiment
                    ev(cc).time = D.time(transTable.onset(trl));
                    ev(cc).duration = 0.05;
                    ev(cc).offset = 0;
                end
            end
            
            % -- GET IMAGES --
            tmp = vertcat(event.stimuli.img{:}); % first get rid of NaNs ('unique' won't work otherwise)
            tmpidx = zeros(length(tmp),1);
            for i = 1:length(tmp)
                if ~isnan(tmp{i,1})
                    tmpidx(i,1) = 1;
                end
            end
            tmp(~tmpidx,:) = [];
            imageNames = unique(tmp);
            imageNames(strcmp(imageNames,'SUPPLY ROOM'),:) = [];
            
            for trl = 1:size(triggerTable,1)
                if any(strcmp(imageNames,triggerTable.label{trl}))
                    
                    % get position of image in sequence
                    stNum = find(strcmp(event.stimuli.img{triggerTable.trial(trl)},...
                                 triggerTable.label{trl}));
                    pathNum = event.stimuli.path{triggerTable.trial(trl)}(stNum,:) + 1; % 1 or 2
                    
                    if ~isempty(stNum) && ~isempty(pathNum)
                        % save to event
                        cc = cc + 1;
                        ev(cc).type = ['image-p' num2str(pathNum) 's' num2str(stNum)];
                        ev(cc).value = triggerTable.trial(trl); % trial number in experiment
                        ev(cc).time = D.time(triggerTable.onset(trl));
                        ev(cc).duration = 0.05;
                        ev(cc).offset = 0;      
                    end
                end
            end
            
            % -- GET OUTCOME -- 
            outIdx = find(strcmp(triggerTable.label,'Outcome') | strcmp(triggerTable.label,'Supply room'));
            outTable = triggerTable(outIdx,:);
            for trl = 1:size(outTable,1)
                cc = cc + 1;
                ev(cc).type = 'outcome';
                ev(cc).value = outTable.trial(trl); % trial number in experiment
                ev(cc).time = D.time(outTable.onset(trl));
                ev(cc).duration = 0.05;
                ev(cc).offset = 0;
            end

        end

        [~,sortidx] = sort(extractfield(ev,'time'));
        ev = ev(sortidx);
        
        disp(struct2table(ev))
        
        D = events(D,1,ev);
        D.save;

        %% Filter

        % files
        opt1 = [];
        opt1.datatype = 'ctf';
        opt1.spm_files = {fullfile(dir_output,['p' subject '_' session '_' num2str(f) '.mat'])};
        opt1.dirname = fullfile(dir_output,'opt1');

        % filter
        opt1.highpass.do = 1;
        opt1.highpass.cutoff = 0.5;
        opt1.mains.do = 1;

        % switch off
        opt1.downsample.do = 0;
        opt1.bad_segments.do = 0;
        opt1.africa.do = 0;
        opt1.outliers.do = 0;
        opt1.coreg.do = 0;
        opt1.epoch.do = 0;

        % run
        opt1 = osl_run_opt(opt1); % saves file starting with 'fffp'

        %% Downsample & bad segments

        % files
        opt2 = [];
        opt2.spm_files = opt1.results.spm_files;
        opt2.datatype = 'ctf';
        opt2.dirname = fullfile(dir_output,'opt2');

        % downsampling
        opt2.downsample.do = 1;
        opt2.downsample.freq = Fs;

        % bad segments
        opt2.bad_segments.do = 1;

        % switch off
        opt2.africa.do = 0;
        opt2.epoch.do = 0;
        opt2.outliers.do = 0;
        opt2.coreg.do = 0;
        opt2.highpass.do = 0;
        opt2.mains.do = 0;

        % run
        opt2 = osl_run_opt(opt2); % saves file starting with 'Bd'
        
        % fix slashes in filename (mac to windows)
        tmp = opt2.results.spm_files{1};
        tmp(strfind(tmp,'/')) = '\';
        opt2.results.spm_files{1} = tmp;

        %% Clean up

        copyfile([opt2.results.spm_files{1} '.mat'],...
                 fullfile(dir_output,['Bdffp_' subject '_' session '_' num2str(f) '.mat']));
        copyfile([opt2.results.spm_files{1} '.dat'],...
                 fullfile(dir_output,['Bdffp_' subject '_' session '_' num2str(f) '.dat']));

        rmdir(fullfile(dir_output,'opt1.opt'),'s');
        rmdir(fullfile(dir_output,'opt2.opt'),'s');
        fnames = dir(fullfile(dir_output));
        for i = 1:length(fnames)
            if fnames(i).name(1) ~= 'B' && fnames(i).name(1) ~= '.'
                delete(fullfile(dir_output,fnames(i).name));
            end
        end
        
    end
    toc
    
end
end
