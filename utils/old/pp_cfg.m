function [filelist,event,in] = pp_cfg(subject,session)

% Get filenames for each subject/sessionand fix any recording errors

% subject  = character (e.g., '088674')
% session  = 'FL' or 'Task'

if isnumeric(subject)
    subject = num2str(subject);
    if length(subject) == 5
        subject = ['0' subject];
    end
    if length(subject) ~= 6
        error(['Invalid subject ID: ' subject])
    end
end

%%

% empty variables
filelist = [];
event = [];
in = [];

% directories
[sinfo, dinfo] = dir_cfg();
dir_behav = fullfile(dinfo.data_behav,subject);
dir_meg = fullfile(dinfo.data_meg_raw,subject);
dir_processed = fullfile(dinfo.data_meg_pp3,subject);

filelist.subject = subject;
filelist.session = session;
filelist.dir_meg = dir_meg;
filelist.dir_processed = dir_processed;

% meg info
in.triggerChannel = 'UADC004';
in.respChannel = 'UPPT002';
in.eog = {'UADC001','UADC003','UADC005'};

% experiment info
sessions = {'FL','Task'}; % possible values for 'session' argument
sessIdx = find(strcmp(sessions,session)); % get index of the current session requested

in.pathStateCodes = [0,0;0,1;0,2;1,0;1,1;1,2]; % path, state
in.negatorCombos = sinfo.nC;

%% subject-specific

switch subject
    case '088674'
        fileID_pre = 'MG06139_RiskyRe_20200124_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '263098'
        fileID_pre = 'MG06142_RiskyRe_20200130_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '680913'
        fileID_pre = 'MG06144_RiskyRe_20200130_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '542599'
        fileID_pre = 'MG06156_RiskyRe_20200204_';
        fileID_post = {
            [1:4]   % FL
            [5:12]  % Task
            };
    case '383991'
        fileID_pre = 'MG06159_RiskyRe_20200205_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '707132'
        fileID_pre = 'MG06164_ClearHD_20200207_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '396430'
        fileID_pre = 'MG06162_RiskyRe_20200206_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '521846'
        fileID_pre = 'MG06175_RiskyRe_20200210_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '015397'
        fileID_pre = 'MG06176_RiskyRe_20200210_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '663186'
        fileID_pre = 'MG06180_RiskyRe_20200212_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '097403'
        fileID_pre = 'MG06199_RiskyRe_20200221_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '503608'
        fileID_pre = 'MG06213_RiskyRe_20200226_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '147947'
        fileID_pre = 'MG06215_RiskyRe_20200227_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    case '506559'
        fileID_pre = 'MG06229_RiskyRe_20200303_';
        fileID_post = {
            [1:4]   % FL
            [5:15]  % Task
            };
    otherwise
        error('Incorrect or missing subject label')
end

behavFiles = {[subject '_fl.csv'],[subject '_task.csv']}; % FL, Task
if ~exist( fullfile(dir_behav,behavFiles{sessIdx}) )
    error([session ' behavioural file is missing for ' subject]);
end

%% session-specific
switch session
    case 'FL' % Functional localiser
        event.labels = {'Image','Cue','Feedback'};
        event.dur = [500,100,300];
    case 'Task'
        event.labels = {'Image','Transition','Choice'};
        event.dur = [500,300,100];
        
        fid = fopen(fullfile(dir_behav,[num2str(subject) '_s2_data.json']));
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        behav = jsondecode(str);
        
        pracStim = behav.x5_negator_learning.stimuli;
        testStim = behav.x6_test.stimuli;
        
        event.stimuli.img = [pracStim.img; testStim.img];
        
        if strcmp(subject,'263098') && strcmp(session,'Task') % for this subject, the practice data wasn't read in as individual cells for some reason
            event.stimuli.state = {};
            event.stimuli.path = {};
            event.stimuli.value = {};
            for i = 1:size(pracStim.state,1)
                event.stimuli.state{i,1} = pracStim.state(i,:)';
                event.stimuli.path{i,1} = pracStim.path(i,:)';
                event.stimuli.value{i,1} = pracStim.value(i,:)';
            end
            for i = 1:size(testStim.state,1)
                event.stimuli.state{i+size(pracStim.state,1),1} = testStim.state{i,1};
                event.stimuli.path{i+size(pracStim.state,1),1} = testStim.path{i,1};
                event.stimuli.value{i+size(pracStim.state,1),1} = testStim.value{i,1};
            end
        else
            event.stimuli.state = [pracStim.state; testStim.state];
            event.stimuli.path = [pracStim.path; testStim.path];
            event.stimuli.value = [pracStim.value; testStim.value];
        end
        
    otherwise
        error('Incorrect or missing session label')
end

%% retrieve data

% get session (FL = 1, Task = 2)
sessIdx = find(not(cellfun('isempty',strfind(sessions,session))));

% get raw MEG files
filelist.data = {};
for i = 1:length(fileID_post{sessIdx})
    filelist.data{i} = [fileID_pre, num2str(fileID_post{sessIdx}(i),'%02.f'), '.ds'];
end

% get behavioural data
event.table = readtable(fullfile(dir_behav,behavFiles{sessIdx}) ,'Delimiter',',','HeaderLines',0);
event.table.Include = ones(size(event.table,1),1);

switch session
    case 'FL' % Functional localiser
        
        % Get block index for each file
        filelist.block_idx = unique(event.table.Block);
        
        % Put in trial numbers
        event.table.Trial = nan(size(event.table,1),1);
        for b = 1:length(filelist.block_idx)
            event.table.Trial(event.table.Block == filelist.block_idx(b)) = 1:sum(event.table.Block == filelist.block_idx(b));
        end
        
        % Get stimulus ID numbers
        stimID = nan(size(event.table,1),1);
        cc = 0;
        for path = 0:1
            for state = 0:2
                idx = event.table.Path == path & event.table.State == state;
                if sum(idx) > 0
                    cc = cc + 1;
                    stimLabel = unique(event.table.Image(event.table.Path == path & event.table.State == state));
                    if length(stimLabel) > 1
                        error('Too many unique image names found!')
                    else
                        stimLabel = stimLabel{1};
                    end
                    stimID(idx) = cc;
                end
            end
        end
        event.table.ID = stimID;
        
    case 'Task'
        
        % Get index of which file is practice vs. test
        block_idx = [event.table.Block, event.table.Practice];
        filelist.block_idx = block_idx(1,:);
        for i = 2:size(block_idx,1)
            if sum(block_idx(i,:) ~= block_idx(i-1,:)) ~= 0
                filelist.block_idx = [filelist.block_idx; block_idx(i,:)];
            end
        end
        
        % Missing behavioural data for subject 680913
        if strcmp(subject,'680913')
            tmp1 = event.table(15,:);
            tmp1(1,:) = array2table(NaN);
            tmp1.Trial = 4;
            tmp1.Choice = 1;
            tmp1.RT = 8.8367;
            tmp1.Block = event.table.Block(15);
            tmp1.Practice = event.table.Practice(15);
            tmp1.ExpTrial = 16;

            tmp2 = event.table(31,:);
            tmp2(1,:) = array2table(NaN);
            tmp2.RT = 5.8967;
            tmp2.Trial = 3;
            tmp2.Choice = 1;
            tmp2.Block = event.table.Block(31);
            tmp2.Practice = event.table.Practice(31);
            tmp2.ExpTrial = 33;
            
            tmp3 = event.table(66,:);
            tmp3(1,:) = array2table(NaN);
            tmp3.RT = NaN;
            tmp3.Trial = 3;
            tmp3.Choice = 1;
            tmp3.Block = event.table.Block(66);
            tmp3.Practice = event.table.Practice(66);
            tmp3.ExpTrial = 69;
            
            tmp4 = event.table(101,:);
            tmp4(1,:) = array2table(NaN);
            tmp4.RT = NaN;
            tmp4.Trial = 3;
            tmp4.Choice = 1;
            tmp4.Block = event.table.Block(101);
            tmp4.Practice = event.table.Practice(101);
            tmp4.ExpTrial = 105;
            
            tmp5 = event.table(103,:);
            tmp5(1,:) = array2table(NaN);
            tmp5.RT = NaN;
            tmp5.Trial = 6;
            tmp5.Block = event.table.Block(103);
            tmp5.Practice = event.table.Practice(103);
            tmp5.ExpTrial = 108;
            
            tmp6 = event.table(118,:);
            tmp6(1,:) = array2table(NaN);
            tmp6.RT = NaN;
            tmp6.Trial = 4;
            tmp6.Choice = 1;
            tmp6.Block = event.table.Block(118);
            tmp6.Practice = event.table.Practice(118);
            tmp6.ExpTrial = 124;
            
            tmp7 = event.table(134,:);
            tmp7(1,:) = array2table(NaN);
            tmp7.RT = NaN;
            tmp7.Trial = 3;
            tmp7.Choice = 1;
            tmp7.Block = event.table.Block(134);
            tmp7.Practice = event.table.Practice(134);
            tmp7.ExpTrial = 141;
            
            tmp8 = event.table(169,:);
            tmp8(1,:) = array2table(NaN);
            tmp8.RT = NaN;
            tmp8.Trial = 3;
            tmp8.Choice = 1;
            tmp8.Block = event.table.Block(169);
            tmp8.Practice = event.table.Practice(169);
            tmp8.ExpTrial = 177;
            
            tmp9 = event.table(170,:);
            tmp9(1,:) = array2table(NaN);
            tmp9.RT = NaN;
            tmp9.Trial = 5;
            tmp9.Block = event.table.Block(170);
            tmp9.Practice = event.table.Practice(170);
            tmp9.ExpTrial = 179;

            event.table = [event.table(1:15,:);
                           tmp1;
                           event.table(16:31,:);
                           tmp2;
                           event.table(32:66,:);
                           tmp3;
                           event.table(67:101,:);
                           tmp4;
                           event.table(102:103,:);
                           tmp5;
                           event.table(104:118,:);
                           tmp6;
                           event.table(119:134,:);
                           tmp7;
                           event.table(135:169,:);
                           tmp8;
                           event.table(170,:);
                           tmp9;
                           event.table(171:end,:)];
                     
             idx = event.table.ExpTrial(isnan(event.table.nCombo));
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.img;
             event.stimuli.img = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.state;
             event.stimuli.state = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.path;
             event.stimuli.path = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.value;
             event.stimuli.value = tmp_img;
             
        elseif strcmp(subject,'097403')
            
            tmp1 = event.table(117,:);
            tmp1(1,:) = array2table(NaN);
            tmp1.Trial = 16;
            tmp1.Block = event.table.Block(117);
            tmp1.Practice = event.table.Practice(117);
            tmp1.ExpTrial = 118;
            
            tmp2 = event.table(180,:);
            tmp2(1,:) = array2table(NaN);
            tmp2.Trial = 8;
            tmp2.Block = event.table.Block(180);
            tmp2.Practice = event.table.Practice(180);
            tmp2.ExpTrial = 182;
            
            tmp3 = event.table(184,:);
            tmp3(1,:) = array2table(NaN);
            tmp3.Trial = 13;
            tmp3.Block = event.table.Block(184);
            tmp3.Practice = event.table.Practice(184);
            tmp3.ExpTrial = 187;
            
            event.table = [event.table(1:117,:);
                           tmp1;
                           event.table(118:180,:);
                           tmp2;
                           event.table(181:184,:);
                           tmp3;
                           event.table(185:end,:)];
                     
             idx = event.table.ExpTrial(isnan(event.table.nCombo));
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.img;
             event.stimuli.img = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.state;
             event.stimuli.state = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.path;
             event.stimuli.path = tmp_img;
             
             tmp_img = cell(size(event.table,1),1);
             tmp_img(idx,:) = {NaN};
             tmp_img(setdiff(1:size(tmp_img,1),idx),:) = event.stimuli.value;
             event.stimuli.value = tmp_img;
            
        end
end

end