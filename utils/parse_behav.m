function results = parse_behav(subject,dir_data)

%% Settings

negators = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]; % negator combos

%% Import data

behav = [];

fbehav = fullfile(dir_data,subject,[num2str(str2num(subject)) '_s2_data.json']);
fstruct = fullfile(dir_data,subject,[num2str(str2num(subject)) '_s2_structure.js']);
fsave = fullfile(dir_data,subject,[num2str(str2num(subject)) '_task.csv']);

% Behavioural output
fid = fopen(fbehav); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
behav = jsondecode(str);

% results = mean(behav.x0_functional_localiser.acc);

% Experiment structure
fid = fopen(fstruct); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
if strcmp(str(1),'v')
    str = str(16:end-2);
end
bstruct = jsondecode(str);

%% Get rows and columns

nTrls = length(struct2array(bstruct.Trial));

colnames = fieldnames(bstruct);
binfo = array2table(nan(nTrls,length(colnames)),'VariableNames',colnames);
for col = 1:5
    thisCol = extractfield(bstruct,colnames{col});
    x = struct2array(thisCol{1})';
    if size(x,2) > 1
        x = num2cell(x,2);
    else
        x = array2table(x);
    end
    binfo(:,col) = x; 
end

%% Parse info

rNames = {'Practice','Block','Trial','Catch','Forced','ExpTrial','nCombo','P','nV_1','nV_2',...
          'S1','S2','S3','S4','S5','S6','nS1','nS2','nS3','nS4','nS5','nS6',...
          'EV','Choice','RT','Acc','Transition','Outcome','Timestamp'};
results = array2table(nan(nTrls,length(rNames)),'VariableNames',rNames);
for trl = 1:nTrls

    % get basic information (trial num, nCombo, P, nV_1, nV_2, EV)
    results(trl,1:5) = binfo(trl,1:5);
    results.ExpTrial(trl) = trl;

    nlist = struct2array(bstruct.N)';
    results.nCombo(trl) = find(negators(:,1) == nlist(trl,1) & negators(:,2) == nlist(trl,2));

    plist = struct2array(bstruct.P)';
    results.P(trl) = plist(trl,1);

    nVlist = reshape(extractfield(bstruct.nV,['x' num2str(trl-1)]),2,3);
    results.nV_1(trl) = nVlist(1,end);
    results.nV_2(trl) = nVlist(2,end);

    evlist = struct2array(bstruct.EV)';
    results.EV(trl) = evlist(trl,1);

    % see if is practice or not
    isPractice = results.Practice(trl)==1;

    % get path values (normal and flipped)
    Vlist = reshape(extractfield(bstruct.V,['x' num2str(trl-1)]),2,3);
    results(trl,11:16) = array2table([Vlist(1,:) Vlist(2,:)]);
    results(trl,17:22) = array2table([nVlist(1,:) nVlist(2,:)]);

    % get subject responses
    if isPractice
        tmp = behav.x5_negator_learning;
    else
        tmp = behav.x6_test;
    end

    idx = find(tmp.trial == extractfield(bstruct.Trial,['x' num2str(trl-1)]) & tmp.block == extractfield(bstruct.Block,['x' num2str(trl-1)]));
    if isempty(idx)
        warning(['Subject ' subject ': could not find trial ' num2str(extractfield(bstruct.Trial,['x' num2str(trl-1)])) ' in block ' ...
            num2str(extractfield(bstruct.Block,['x' num2str(trl-1)]))])
        results(trl,:) = array2table(nan(1,size(results,2)));
    else

        if ~isempty(strmatch(tmp.choice{idx},'Airlock'))
            results.Choice(trl) = 1;
        elseif ~isempty(strmatch(tmp.choice{idx},'SUPPLY ROOM'))
            results.Choice(trl) = 2;
        end

        results.RT(trl) = tmp.rt(idx)/1000; % in seconds
        results.Acc(trl) = (results.Choice(trl) == 1 && results.EV(trl) >= 1) || (results.Choice(trl) == 2 && results.EV(trl) <= 1);
        if extractfield(bstruct.Forced,['x' num2str(trl-1)]) ~= 0
            results.Acc(trl) = NaN;
        end

        if results.Choice(trl) == 2
            results.Outcome(trl) = 1;
        else
            if ~isempty(strmatch(tmp.transition{idx},'DOOR 1'))
                results.Outcome(trl) = results.nV_1(trl);
            elseif ~isempty(strmatch(tmp.transition{idx},'DOOR 2'))
                results.Outcome(trl) = results.nV_2(trl);
            end
        end

        if strcmp(tmp.transition(idx),'DOOR 1')
            results.Transition(trl) = 1;
        elseif strcmp(tmp.transition(idx),'DOOR 2')
            results.Transition(trl) = 2;
        elseif strcmp(tmp.transition(idx),'SUPPLY ROOM')
            results.Transition(trl) = 0;
        end
        
        results.Timestamp(trl) = tmp.timeStamp(idx);
        
    end
end

% remove all NaN rows
idx = [];
for trl = 1:size(results,1)
    tmp = isnan(table2array(results(trl,:)));
    if sum(tmp) == length(tmp)
        idx = [idx; trl];
    end
end
results(idx,:) = [];

% Clean up other factors
results.bAcc = results.Acc; % put 1 point buffer either side of EV to score accuracy
results.bAcc(results.EV >= 0 & results.Choice == 1) = 1;
results.bAcc(results.EV <= 2 & results.Choice == 2) = 1;
results.bAcc(isnan(results.Acc)) = NaN;
results.RT = results.RT + 5; % the first 5 seconds aren't included in the RT

% Manual fix for errors caused by non-responses
if strcmp(subject,'097403')
    results(results.Practice==0 & results.Block==6 & (results.Trial==15 | results.Trial==17),:) = [];
    results(results.Practice==0 & results.Block==10 & (results.Trial==7 | results.Trial==12),:) = [];
end
if strcmp(subject,'680913')
    results(results.Practice==0 & results.Block==1 & (results.Trial==1 | results.Trial==3),:) = [];
    results(results.Practice==0 & results.Block==2 & (results.Trial==2 | results.Trial==3),:) = [];
    results(results.Practice==0 & results.Block==4 & (results.Trial==2),:) = [];
    results(results.Practice==0 & results.Block==6 & (results.Trial==2 | results.Trial==4 | results.Trial==5),:) = [];
    results(results.Practice==0 & results.Block==7 & results.Trial<5,:) = [];
    results(results.Practice==0 & results.Block==8 & (results.Trial==2),:) = [];
    results(results.Practice==0 & results.Block==10 & results.Trial<6,:) = [];
end
if strcmp(subject,'945092')
    results(results.Practice==1 & results.Trial<=4,:) = []; % had to start MEG from trial 5 onwards because participant asked question during trial 4
end

% Write table
results.Subject = repmat(subject,size(results,1),1);

disp(['Subject ' subject ': ' num2str(size(results,1)) ' trials in total ',...
    '(' num2str(sum(~isnan(results.Acc))) ' non-nan), ',...
    num2str(round(nanmean(results.Acc)*100,2)) '% accuracy'])

writetable(results,fsave);

end