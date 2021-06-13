function func_epoch(datafile,params)

tic

[subject,session,event,in] = pp_osl_cfg(datafile);
disp(['----- FILENAME: ' datafile])
disp(['----- SUBJECT:  ' subject])
disp(['----- SESSION:  ' session])

% update output directory
params.dir_output = fullfile('/Volumes/SOFTWARE BA/processed',subject);
addpath(params.dir_output);

%% Epoch

fparts = strsplit(datafile,'_');
dtype = fparts{3}; % 'FL' or 'Task'
fnum = strsplit(fparts{end},'.mat');
fnum = str2num(fnum{1});

ica_file = fullfile(params.dir_output,[subject '_' dtype '_' num2str(fnum) '_ICA.mat']);

opt4 = [];

opt4.spm_files = {ica_file};
opt4.datatype = 'ctf';
opt4.dirname = fullfile(params.dir_output,'opt4');

% epoch settings
opt4.epoch.do = 1;
opt4.outliers.do = 1;

% coreg settings
opt4.coreg.do = 1;
opt4.coreg.use_rhino = 0;
opt4.coreg.useheadshape = 0;

% switch off
opt4.bad_segments.do = 0;
opt4.downsample.do = 0;
opt4.africa.todo.ica = 0;
opt4.africa.todo.ident = 0;
opt4.africa.todo.remove = 0;

% --- Define trials using photodiode output
if strcmp(session,'FL')
    opt4.epoch.time_range = [-.1 .5];
    for st = 1:length(unique(event.table.Image))
        opt4.epoch.trialdef(st).conditionlabel = event.table.Image{find(event.table.ID == st,1,'first')};
        opt4.epoch.trialdef(st).eventtype = 'image';
        opt4.epoch.trialdef(st).eventvalue = st;
    end
elseif strcmp(session,'Task')
    
    D = spm_eeg_load(opt4.spm_files{1});
    ev = D.events;
    
    if strcmp(params.epoch,'decision')

        all_vals = nan(0,2);
        cc = 0;
        for e = 1:length(ev)
            if strcmp(ev(e).type,'choice')
                cc = cc + 1;
                all_vals(cc,1) = ev(e).value;
                all_vals(cc,2) = event.table.Choice(event.table.ExpTrial == ev(e).value);
            end
        end

        opt4.epoch.trialdef(1).conditionlabel = 'risky';
        opt4.epoch.trialdef(1).eventtype = 'choice';
        opt4.epoch.trialdef(1).eventvalue = all_vals(all_vals(:,2) == 1,1);

        opt4.epoch.trialdef(2).conditionlabel = 'safe';
        opt4.epoch.trialdef(2).eventtype = 'choice';
        opt4.epoch.trialdef(2).eventvalue = all_vals(all_vals(:,2) == 2,1);
        
        maxRT = ceil(max(event.table.RT(event.table.Block == in.block & event.table.Practice == in.practice)));

        opt4.epoch.time_range = [-.5 maxRT]; % take maximum decision window (to trim later to RTs)
        
    elseif strcmp(params.epoch,'transition')
        
        opt4.epoch.time_range = [-.5 3];
        
        all_vals = nan(0,2);
        cc = 0;
        for e = 1:length(ev)
            if strcmp(ev(e).type,'transition')
                cc = cc + 1;
                all_vals(cc,1) = ev(e).value;
                all_vals(cc,2) = event.table.Transition(event.table.ExpTrial == ev(e).value);
            end
        end

        opt4.epoch.trialdef(1).conditionlabel = 'transition1';
        opt4.epoch.trialdef(1).eventtype = 'transition';
        opt4.epoch.trialdef(1).eventvalue = all_vals(all_vals(:,2) == 1,1);

        opt4.epoch.trialdef(2).conditionlabel = 'transition2';
        opt4.epoch.trialdef(2).eventtype = 'transition';
        opt4.epoch.trialdef(2).eventvalue = all_vals(all_vals(:,2) == 2,1);
        
    elseif strcmp(params.epoch,'image')
        
        % check image labels against event log
        allvals = cell(2,3);
        for i = 1:length(ev)
            if ~isempty(strfind(ev(i).type,'image'))
                thisPath = str2num(ev(i).type(8));
                thisState = str2num(ev(i).type(10));
                if event.table.Transition(event.table.ExpTrial == ev(i).value) ~= thisPath
                    error('Paths do not match')
                else
                   allvals{thisPath,thisState} = [allvals{thisPath,thisState}; ev(i).value];
                end
            end
        end
        allvals = reshape(allvals',1,6);

        opt4.epoch.time_range = [-.2 1];
       
        % get stimuli names
        stimNames = cell(1,6);
        for i = 1:length(event.stimuli.img)
            if unique(event.stimuli.path{i}) == 0 && sum(cellfun('isempty',stimNames(1:3))) == 3
                stimNames(1:3) = event.stimuli.img{i}';
            elseif unique(event.stimuli.path{i}) == 1 && sum(cellfun('isempty',stimNames(4:6))) == 3
                stimNames(4:6) = event.stimuli.img{i}';
            end
            if ~any(cellfun('isempty',stimNames))
                break
            end
        end
        
        for i = 1:6
            opt4.epoch.trialdef(i).conditionlabel = stimNames{i};
            if i <= 3
                opt4.epoch.trialdef(i).eventtype = ['image-p1s' num2str(i)];
            elseif i >= 4
                opt4.epoch.trialdef(i).eventtype = ['image-p2s' num2str(i-3)];
            end
            opt4.epoch.trialdef(i).eventvalue = allvals{i};
        end
        
    elseif strcmp(params.epoch,'outcome')
 
        opt4.epoch.time_range = [-.5 2.5];
        
        all_vals = nan(0,3);
        cc = 0;
        for e = 1:length(ev)
            if strcmp(ev(e).type,'outcome')
                cc = cc + 1;
                all_vals(cc,1) = ev(e).value;
                all_vals(cc,2) = event.table.Outcome(event.table.ExpTrial == ev(e).value);
                all_vals(cc,3) = event.table.Choice(event.table.ExpTrial == ev(e).value);
            end
        end
        
        opt4.epoch.trialdef(1).conditionlabel = 'supplyroom';
        opt4.epoch.trialdef(1).eventtype = 'outcome';
        opt4.epoch.trialdef(1).eventvalue = all_vals(all_vals(:,3) == 2,1);

        uVals = unique(all_vals(all_vals(:,3) == 1,2)); % get unique outcomes for risk choices
        for u = 1:length(uVals)
            if uVals(u) > 0
                opt4.epoch.trialdef(u+1).conditionlabel = ['+' num2str(uVals(u))];
            elseif uVals(u) < 0
                opt4.epoch.trialdef(u+1).conditionlabel = ['' num2str(uVals(u))];
            else
                error('Outcome value was zero?')
            end
            opt4.epoch.trialdef(u+1).eventtype = 'outcome';
            opt4.epoch.trialdef(u+1).eventvalue = all_vals(all_vals(:,2) == uVals(u),1);
        end
        
    end    

end

% run
opt4 = osl_run_opt(opt4);

if isempty(opt4.results.spm_files_epoched{1})
    error(['Epochs did not run properly for file ' datafile])
end

%% Clean up

finalFile = fullfile(params.dir_output,[subject '_' session '_' num2str(in.idx) '_' params.epoch '.mat']);

D = spm_eeg_load(opt4.results.spm_files_epoched{1});
dinfo = [];
dinfo.badtrials = D.badtrials;
dinfo.badchannels = D.badchannels;
dinfo.epochinfo = D.epochinfo;
dinfo.events = D.events;
dinfo.fid = D.fiducials;
dinfo.ica = D.ica;

D = spm2fieldtrip(D);
save(finalFile,'D');
save(fullfile(params.dir_output,[subject '_' session '_' num2str(in.idx) '_' params.epoch '_info.mat']),'dinfo');

try
    rmdir(fullfile(params.dir_output,'opt3.opt'),'s');
end
try
    rmdir(fullfile(params.dir_output,'opt4.opt'),'s');
end

copyfile(fullfile('/Volumes/SOFTWARE BA/halfprocessed/todo/',[datafile(1:end-4) '.mat']),...
    fullfile('/Volumes/SOFTWARE BA/halfprocessed/done/',[datafile(1:end-4) '.mat']));
copyfile(fullfile('/Volumes/SOFTWARE BA/halfprocessed/todo/',[datafile(1:end-4) '.dat']),...
    fullfile('/Volumes/SOFTWARE BA/halfprocessed/done/',[datafile(1:end-4) '.dat']));
delete(fullfile('/Volumes/SOFTWARE BA/halfprocessed/todo/',[datafile(1:end-4) '.mat']));
delete(fullfile('/Volumes/SOFTWARE BA/halfprocessed/todo/',[datafile(1:end-4) '.dat']));

toc

end