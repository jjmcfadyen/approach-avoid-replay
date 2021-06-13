% Main preprocessing script

clear all
clc

%% Directories

addpath('D:\Toolboxes\spm12')
spm('defaults','eeg')

dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';
dir_output = 'D:\2020_RiskyReplay\data\meg';
dir_behav = 'D:\2020_RiskyReplay\data\behav';

%% Parameters

addpath('utils');
addpath('preprocessing')

parameters = get_parameters(dir_raw);

% remove practice blocks for 018768 & 957849 (photodiode missing)
parameters((strcmp(parameters.schar,'018768') | strcmp(parameters.schar,'957849')) & ...
    strcmp(parameters.task,'task') & parameters.block==0,:) = [];

% remove block 9 for 391883 (photodiode signal lost)
parameters(strcmp(parameters.schar,'391883') & strcmp(parameters.task,'task') & parameters.block==9,:) = [];

subjects = unique(parameters.schar);
N = length(subjects);

%% Convert data to SPM

for s = 1:N
    
    % Get list of raw data files
    cd(fullfile(dir_raw,subjects{s}))
    idx = find(parameters.subjectID==str2double(subjects{s}));
    filelist = parameters.rawfile(idx);
    
    % Set output directory
    thisoutput = fullfile(dir_output,'1_converted',subjects{s});
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
    
    % Set input directory
    thisinput = fullfile(dir_output,'1_converted',subjects{s});
    if ~exist(thisinput)
        mkdir(thisinput)
    end
    
    % Set output directory
    thisoutput = fullfile(dir_output,'2_cropped',subjects{s});
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

figure
set(gcf,'position',[5 219 1431 777])
for s = 1:N
   
    dir_input = fullfile(dir_output,'2_cropped',subjects{s});
    
    D = spm_eeg_load(fullfile(dir_input,...
        ['cropped_' subjects{s} '_' thistask '_r1.mat']));
    
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

    D = chantype(D,find(contains(D.chanlabels,eog)),'EOG');
    D.save;
    
end

%% Generate jobs for OSL on cluster (filter/ICA/epoch)

addpath('D:\2020_RiskyReplay\approach-avoid-replay\preprocessing\batch')

for s = 1:N
   
    filelist = dir(fullfile(dir_output,'2_cropped',subjects{s},'cropped*.mat'));
    
    for f = 1:length(filelist)
        generate_jobs_preprocess(filelist(f).name);
    end
    
end

