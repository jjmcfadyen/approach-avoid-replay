function [triggers,endtime] = parse_triggers(D,behav,task)

subject = unique(behav.Subject);

%% Settings

pLabel = 'UADC004'; % channel label for photodiode

if strcmp(subject,'018768') || strcmp(subject,'957849')
    pLabel = 'UADC005';
end

%% Extract photodiode

pData = squeeze(D(find(contains(D.chanlabels,pLabel)),:,:));

%% Get characteristics of data trace

croptime = find(pData==0,1,'first')-1;
disp(['Photodiode trace ends after ' num2str(round(D.time(croptime)/60,2)) ' minutes'])

if ~strcmp(subject,'957849')
    pData = detrend(pData(1:croptime),1);
else
    pData = pData(1:croptime);
end
baseline = median(pData); % most common value is the baseline amplitude
thisdev = std(pData)/2; % with some fluctation

if strcmp(subject,'663186') && strcmp(task,'task') && all(behav.Block==7)
   thisdev =  std(pData)*1.5; % screen kicked, some blips arehigher amplitude
elseif strcmp(subject,'957849') && strcmp(task,'FL') && all(behav.Block==1)
   thisdev =  .001; % amplitude dims in last third of block?
end

onsets = pData > baseline+thisdev;

% Get onsets/offsets/durations
triggers = [];
thistrig = [];
for i = 2:length(onsets)
   
    if onsets(i)==1 && onsets(i-1)==0 % 0 to 1
        thistrig.sOnset = i;
        thistrig.tOnset = D.time(i);
    elseif onsets(i)==0 && onsets(i-1)==1 % 1 to 0
        thistrig.sOffset = i;
        thistrig.tOffset = D.time(i);
        thistrig.sDur = thistrig.sOffset-thistrig.sOnset;
        thistrig.tDur = thistrig.tOffset-thistrig.tOnset;
        triggers = [triggers; struct2table(thistrig)];
    end
    
end

% Clean up any obvious errors
triggers(triggers.tDur < .04,:) = []; % remove any blips shorter than 40 ms
triggers(triggers.tDur > .7,:) = [];  % remove any blips longer than 700 ms

% Make manual fixes
if strcmp(subject,'220598') && strcmp(task,'FL') && all(behav.Block==2)
    triggers.tDur(77) = .1;
end
if strcmp(subject,'088674') && strcmp(task,'task') && all(behav.Block==9)
    triggers(2,:) = [];
end
if strcmp(subject,'263098') && strcmp(task,'task') && all(behav.Practice==0)
    if all(behav.Block==9)
        triggers(2:3,:) = [];
    elseif all(behav.Block==10)
        triggers(2,:) = [];
    end
end
if strcmp(subject,'503608') && strcmp(task,'task') && all(behav.Practice==0)
    if all(behav.Block==1)
        triggers(triggers.tOnset > 203.7 & triggers.tOnset < 208,:) = [];
    end
end
if strcmp(subject,'957849') && strcmp(task,'FL') && all(behav.Block==3)
    triggers(end+1,:) = array2table([20089 16.74 20677 17.23 501 .49],'variablenames',triggers.Properties.VariableNames);
    triggers(end+1,:) = array2table([21301 17.75 21421 17.85 102 .1],'variablenames',triggers.Properties.VariableNames);
    triggers(end+1,:) = array2table([22465 18.72 22813 19.01 296 .29],'variablenames',triggers.Properties.VariableNames);
    triggers = sortrows(triggers);
end

%% Match to behavioural log

switch task
    case 'FL'
        
        % --------------------------
        % Triggers:
        % --------------------------
        % 1. Picture onset  (500 ms)
        % 2. Word prompt    (100 ms)
        % 3. Fixation cross (300 ms) - these are ommitted because they're not present for all subjects
        % --------------------------
        
        triggerNames = {'picture','words'};
        triggerDur = [.5 .1];
        triggers(round(triggers.tDur,1)==.3,:) = []; % ignoring fixation cross triggers as they're missing for many participants
        
        nTrigs = size(triggers,1);
        nTypes = length(triggerNames);
        
        % Check overall trigger number matches (there should be 3 triggers per trial)
        nTrls = size(behav,1);
        
        if size(triggers,1) ~= nTypes*nTrls
            error('Trigger mismatch')
        end
        
        % Then see that the order of triggers makes sense (it should go 1,2,3 for each trial)
        triggerOrder = round(triggers.tDur,1); % trigger type, according to duration
        for i = 1:nTypes
            triggerOrder(triggerOrder==triggerDur(i),2) = i;
        end
        
        if any(triggerOrder(:,2) ~= repmat([1:nTypes]',nTrls,1))
            error('Trigger mismatch')
        end
        
        % If no errors, then add descriptors
        triggers.type = cell(nTrigs,1);
        triggers.trial = nan(nTrigs,1);
        triggers.block = nan(nTrigs,1);
        triggers.acc = nan(nTrigs,1);
        triggers.rt = nan(nTrigs,1);
        triggers.path = nan(nTrigs,1);
        triggers.state = nan(nTrigs,1);
        triggers.imgName = cell(nTrigs,1);
        
        for i = 1:nTypes
            triggers.type(triggerOrder(:,2)==i) = triggerNames(i);
            triggers.trial(i:nTypes:end) = behav.Trial;
            triggers.block(i:nTypes:end) = behav.Block;
            triggers.acc(i:nTypes:end) = strcmp(behav.Acc,'True');
            triggers.rt(i:nTypes:end) = behav.RT;
            triggers.path(i:nTypes:end) = behav.Path;
            triggers.state(i:nTypes:end) = behav.State;
            triggers.imgName(i:nTypes:end) = behav.Image;
        end
        
        triggers.imgNum = nan(nTrigs,1);
        cc = 0;
        for p = 0:1 % for each path
            for st = 0:2 % for each state 
                cc = cc + 1;
                triggers.imgNum(triggers.path==p & triggers.state==st) = cc;
            end
        end
        
    case 'task'
        
        % -----------------------------------------
        % Triggers:
        % -----------------------------------------
        % 1. Trial start                  (100 ms)
        % 2. Transition screen            (300 ms)
        % 3. Image onset / final outcome  (500 ms)
        % -----------------------------------------
        
        triggerNames = {'decision','transition','image','outcome'};
        
        nTrigs = size(triggers,1);
        nTypes = length(triggerNames);
        
        % Some trials missing from behavioural log (but present in MEG)
        if strcmp(subject,'097403') && strcmp(task,'task')
            if all(behav.Practice==0) && all(behav.Block==6)
                triggers(73:81,:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Practice==0) && all(behav.Block==10)
                triggers([33:35 54:56],:) = [];
                nTrigs = size(triggers,1);
            end
        end
        if strcmp(subject,'680913') && strcmp(task,'task') && all(behav.Practice==0)
            if all(behav.Block==1)
                triggers([1:5 12:13],:) = [];
                triggers(7:8,:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==2)
                triggers(7:9,:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==4)
                triggers([7:9 76:77],:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==6)
                triggers(7:14,:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==7)
                triggers([1:12 27],:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==8)
                triggers([7:9 30:31],:) = [];
                nTrigs = size(triggers,1);
            elseif all(behav.Block==10)
                triggers(1:13,:) = [];
                nTrigs = size(triggers,1);
            end
        end
        
        % Work out how many triggers there should be, given behavioural log
        startTrigs = size(behav,1); % this many trial starts
        approachTrigs = sum(behav.Choice==1)*5; % 1 x transition, 3 x images, 1 x outcome
        avoidTrigs = sum(behav.Choice==2); % 1 x onset
        
        triggerOrder = round(triggers.tDur,1); % trigger type, according to duration
        triggerOrder(triggerOrder==.1,2) = 1; % trial start
        triggerOrder(triggerOrder==.3,2) = 2; % transition
        triggerOrder(triggerOrder==.5,2) = 3; % image / outcome
        
        behavOrder = [];
        for trl = 1:size(behav,1) % for each trial in behavioural log...
             if behav.Choice(trl)==1 % approach
                 behavOrder = [behavOrder; 1 trl; 2 trl; 3 trl; 3 trl; 3 trl; 3 trl];
             elseif behav.Choice(trl)==2 % avoid
                 behavOrder = [behavOrder; 1 trl; 3 trl];
             end
        end
        
        behavSum = sum([startTrigs approachTrigs avoidTrigs]);
        if nTrigs ~= behavSum
            
            % sometimes there is noise in the photodiode at the end
            if nTrigs > behavSum
                if all(triggerOrder(1:behavSum,2) == behavOrder(:,1))
                    triggers = triggers(1:behavSum,:);
                    triggerOrder = triggerOrder(1:behavSum,:);
                    nTrigs = size(triggers,1);
                else
                    error('Trigger mismatch')
                end
            else
                error('Trigger mismatch')
            end
            
            %{
            figure; plot(D.time(1:croptime),pData); hold on
            scatter(triggers.tOnset,pData(triggers.sOnset));
            round(triggers.tDur,1)*10
            %}
        end
        
        if any(triggerOrder(:,2) ~= behavOrder(:,1))
            error('Trigger mismatch')
        end
        
        % If no errors, mark the 'outcome' screen (has same duration as images)
        triggerOrder(find(diff(triggerOrder(:,2))==-2),2) = 4;
        
        % Add descriptors
        triggers.type = cell(nTrigs,1);
        for i = 1:nTypes
            triggers.type(triggerOrder(:,2)==i) = triggerNames(i);
        end

        for trl = 1:size(behav,1)
            triggers.practice(behavOrder(:,2)==trl) = behav.Practice(trl); 
            triggers.block(behavOrder(:,2)==trl) = behav.Block(trl); 
            triggers.trial(behavOrder(:,2)==trl) = behav.Trial(trl); 
            triggers.forced(behavOrder(:,2)==trl) = behav.Forced(trl)~=0; 
            triggers.ev(behavOrder(:,2)==trl) = behav.EV(trl); 
            triggers.choice(behavOrder(:,2)==trl) = behav.Choice(trl); 
            triggers.rt(behavOrder(:,2)==trl) = behav.RT(trl); 
            triggers.outcome(behavOrder(:,2)==trl) = behav.Outcome(trl); 
            triggers.transition(behavOrder(:,2)==trl) = behav.Transition(trl); 
        end
        
        tmp = cell(nTrigs,1);
        tmp(triggers.choice==1) = {'approach'};
        tmp(triggers.choice==2) = {'avoid'};
        triggers.choice = tmp;
       
        
end

%% Establish end point for cropping

endtime = triggers.tOffset(end) + 5; % final trigger, plus 5 seconds
if endtime > D.time(croptime)
    endtime = D.time(croptime);
    warning(['Recording finishes only ' num2str(D.time(croptime)-triggers.tOffset(end)) ' seconds after last trigger'])
end

end