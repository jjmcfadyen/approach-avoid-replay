function [triggers, fsave] = pp_triggers(fl,event,in,f,bC)

% Finds triggers and matches them to behavioural log
%
%       subject    = character (e.g., '088674')
%       session    = 'FL' or 'Task'
%       bC         = variable produced by pp_findEndSamples
%
% (requires fieldtrip to be added to path)

%% Read each block separately

subject = fl.subject;
session = fl.session;

fprintf(['---------------------------------------------------\n'])
fprintf(['---------------------------------------------------\n'])
fprintf(['------------------ ' subject ' -------------------------\n'])
fprintf(['-------------- (file ' num2str(f) ' of ' num2str(length(fl.data)) ') ---------------------\n'])
fprintf(['---------------------------------------------------\n'])
fprintf(['---------------------------------------------------\n'])

%% Import data

% load raw meg data file
cfg = [];
cfg.dataset = fullfile(fl.dir_meg,fl.data{f});
data = ft_preprocessing(cfg);

Fs = data.fsample;

% crop to end point
data.time{1} = data.time{1}(:,1:bC);
data.trial{1} = data.trial{1}(:,1:bC);
data.sampleinfo = [1 bC];

% select trigger channel
cfg = [];
cfg.channel = in.triggerChannel;
trigger = ft_selectdata(cfg,data);

%     cfg = [];
%     cfg.detrend = 'yes';
%     trigger = ft_preprocessing(cfg,trigger); % remove linear trend & add filter (helps with calculating median value later)

trigger = trigger.trial{1}; % just take the data (ignore other meta info)

% select response button channel
cfg = [];
cfg.channel = in.respChannel;
response = ft_selectdata(cfg,data);
response = response.trial{1}; % just take the data (ignore other meta info)

%% Identify trigger onsets

% mark data that goes above threshold
window = 10; % in seconds
tMedian = movmedian(trigger,window*data.fsample); % get median value
tMax = max(trigger); % get maximum value
tLim = tMedian + (tMax - tMedian) * .1;
spikes = trigger >= tLim;

% get abrupt changes
idx = find(ischange(double(spikes)));

% identify whether the abrupt change is an onset or offset
onsets = [];
offsets = [];
changerate = [NaN diff(trigger)];
for i = 1:length(idx)
    if changerate(idx(i)) > 0
        onsets = [onsets; idx(i)];
    elseif changerate(idx(i)) < 0 
        offsets = [offsets; idx(i)];
    end
end

% move as close as possible to baseline
if strcmp(session,'FL')
    for i = 1:length(onsets)
        onsets(i) = find(trigger(1:onsets(i)) < tMedian(onsets(i)) + abs(tMedian(onsets(i))),1,'last');
    end
    for i = 1:length(offsets)
        offsets(i) = find(trigger(offsets(i):end) < tMedian(offsets(i)) + abs(tMedian(offsets(i))),1,'first') + offsets(i);
    end
elseif strcmp(session,'Task')
    for i = 1:length(onsets)
        onsets(i) = find(trigger(1:onsets(i)) <= tMedian(onsets(i)),1,'last');
    end
    for i = 1:length(offsets)
        offsets(i) = find(trigger(offsets(i):end) <= tMedian(offsets(i)),1,'first') + offsets(i);
    end
end

% delete duplicates
onsets = unique(onsets);
offsets = unique(offsets);

if strcmp(subject,'506559') && strcmp(session,'FL') && f == 3
    offsets(239,:) = [];
elseif strcmp(subject,'680913') && strcmp(session,'Task')
    if f == 9
        onsets([8 31:33],:) = [];
        offsets([8 31 32],:) = [];
    end
elseif strcmp(subject,'503608') && strcmp(session,'Task')
    if f == 2
        offsets(58,:) = []; 
    end
end

if length(onsets) ~= length(offsets)
    error('')
end

% delete impossibly fast triggers
dur = (offsets - onsets)/1000;
idx = dur < (1/60);
onsets(idx,:) = [];
offsets(idx,:) = [];

% manual error corrections
if strcmp(subject,'088674') && strcmp(session,'Task')
    if f == 3
        offsets(9) = 5.67e4;
    elseif f == 4
        offsets(15) = 7.722e4;
        onsets(36) = 1.773e5;
    elseif f == 6
        onsets(41) = 2.036e5;
    elseif f == 10
        onsets(2,:) = [];
        offsets(2,:) = [];
    end
elseif strcmp(subject,'263098') && strcmp(session,'Task')
    if f == 1
         offsets(3) = 2.018e4;
         offsets(46) = 2.733e5;
         onsets(46) = 2.727e5;
    elseif f == 2
        offsets(9) = 4.03e4;
    elseif f == 3
        offsets(21) = 8.778e4;
        onsets(25) = 1.001e5;
        offsets(25) = 1.003e5;
    elseif f == 8
        offsets(9) = 3.871e4;
    elseif f == 10
        offsets(34) = 1.36e5;
    elseif f == 11
        onsets(2,:) = [];
        offsets(2,:) = [];
    end

elseif strcmp(subject,'680913') && strcmp(session,'Task')
    if f == 1
        offsets(15) = 1.015e5;
        offsets(21) = 1.53e5;
    elseif f == 2
        onsets(13:14,:) = []; 
        offsets(13:14,:) = [];
    elseif f == 3
        offsets(13,:) = 9.848e4;
        onsets([6 14],:) = [];
        offsets([6 14],:) = [];
        onsets(13,:) = [];
        offsets(13,:) = [];
    elseif f == 4
        offsets(3,:) = 1.212e4;
        offsets(21,:) = 8.508e4;
    elseif f == 6
        offsets(9,:) = 3.891e4;
    elseif f == 7
        offsets(3,:) = 3.184e4;
    elseif f == 8
        onsets([2 4 28],:) = [];
        offsets([2 4 28],:) = [];
    elseif f == 10
        onsets(16,:) = [];
        offsets(15,:) = 7.403e4;
        offsets(16,:) = [];
    elseif f == 11
        onsets(10,:) = [];
        offsets(10,:) = [];
    end
elseif strcmp(subject,'383991') && strcmp(session,'Task')
    if f == 3
        onsets(14,:) = 8.436e4;
        onsets(31,:) = 1.658e5;
        offsets(31,:) = 1.659e5;
        onsets(43,:) = 2.476e5;
        offsets(43,:) = 2.482e5;
    end
elseif strcmp(subject,'707132') && strcmp(session,'Task')
    if f == 1
        offsets(9,:) = 2.176e4; 
    elseif f == 5
        offsets(15,:) = 4.687e4;
    end
elseif strcmp(subject,'396430') && strcmp(session,'Task')
    if f == 1
        offsets(9,:) = 8.458e4;
    elseif f == 2
        offsets(21,:) = 1.144e5;
        onsets(53,:) = 2.746e5;
        offsets(53,:) = 2.748e5;
    elseif f == 4
        offsets(54,:) = 3.265e5;
    elseif f == 6
        onsets(25,:) = 1.099e5;
        onsets(88,:) = 4.524e5;
        onsets(94,:) = 5.067e5;
        offsets(3,:) = 1.573e4;
        offsets(25,:) = 1.1e5;
        offsets(85,:) = 4.412e5;
        offsets(88,:) = 4.531e5;
        offsets(94,:) = 5.073e5;
    elseif f == 7
        onsets(43,:) = 2.62e5;
        offsets(43,:) = 2.622e5;
    elseif f == 8
        onsets(38,:) = 1.932e5;
        offsets(38,:) = 1.938e5;
    elseif f == 9
        offsets(3,:) = 1.27e4;
    elseif f == 11
        onsets(43,:) = 2.199e5;
        offsets(43,:) = 2.205e5;
    end
elseif strcmp(subject,'521846') && strcmp(session,'Task')
    if f == 3
        offsets(15,:) = 9.311e4; 
    elseif f == 9
        onsets(44,:) = 2.114e5;
    elseif f == 11
        onsets(80,:) = 4.491e5;
    end
elseif strcmp(subject,'015397') && strcmp(session,'Task')
    if f == 2
        onsets(43,:) = 1.746e5;
    elseif f == 3
        onsets(9,:) = 3.827e4;
        offsets(21,:) = 8.558e4;
    end
elseif strcmp(subject,'663186') && strcmp(session,'Task')
    if f == 8
        onsets(28,:) = 1.296e5;
        onsets(29,:) = 1.326e5;
        offsets(28,:) = 1.302e5;
    elseif f == 10
        onsets(49,:) = 2.241e5;
        offsets(49,:) = 2.247e5;
    end
elseif strcmp(subject,'503608') && strcmp(session,'Task')
    if f == 2
        offsets(3,:) = 1.523e4;
        onsets(44,:) = 1.92e5;
        onsets(58,:) = [];
        offsets(58,:) = [];
    elseif f == 3
        offsets(9,:) = 4.269e4;
        onsets(20,:) = 8.917e4;
        offsets(20,:) = 8.954e4;
        offsets(21,:) = 9.28e4;
        onsets(69,:) = 3.661e5;
        onsets(78,:) = 4.117e5;
        onsets(87,:) = 4.554e5;
        offsets(87,:) = 4.56e5;
    elseif f == 4
        offsets(3,:) = 1.592e4;
    elseif f == 7
        offsets(54,:) = 2.558e5;
    end
elseif strcmp(subject,'147947') && strcmp(session,'Task')
    if f == 1
        onsets(49,:) = 2.393e5;
    elseif f == 11
        offsets(15,:) = 6.402e4;
    end
elseif strcmp(subject,'506559') && strcmp(session,'Task')
    if f == 4
        offsets(15,:) = 9.714e4;
    end
end

if length(onsets) ~= length(offsets)
    error([subject ', ' session ' file ' num2str(f) ': different no. of trigger onsets & offsets'])
end

% put into table
triggers = array2table([onsets offsets offsets-onsets],'VariableNames',{'onset','offset','dur'});

% plot
figure
plot(trigger,'k'); hold on
scatter(onsets,trigger(onsets),50,'g'); hold on
scatter(offsets,trigger(offsets),50,'r'); hold on
plot(tMedian,'g'); hold on
plot(tLim,'b'); hold on
title('Photodiode')
xlabel('Samples')
ylabel('Signal amplitude')

%% Identify button presses

changerate = [0 diff(response)];
idx = find(changerate ~= 0); % find where the response signal changes
onsets = [];
offsets = [];
for i = 1:length(idx)
    if changerate(idx(i)) > 0
        onsets = [onsets; [idx(i) abs(changerate(idx(i)))]];
    else
        offsets = [offsets; [idx(i) abs(changerate(idx(i)))]];
    end
end

% put into table
responses = array2table([onsets(:,1) offsets],'VariableNames',{'onset','offset','type'});

%% Label triggers

% use trigger durations to guess labels
dur = triggers.dur / data.fsample;

[~,~,bin] = histcounts(dur,length(event.dur));
if strcmp(subject,'707132') && strcmp(session,'Task') && f == 11
    bin(dur < .15,:) = 1;
    bin(dur > .15 & dur < .5,:) = 2;
    bin(dur > .5,:) = 3;
elseif strcmp(subject,'396430') && strcmp(session,'Task')
    if f == 4
        bin(dur < .15,:) = 1;
        bin(dur > .15 & dur < .5,:) = 2;
        bin(dur > .49,:) = 3;
    elseif f == 5
        bin(dur < .15,:) = 1;
        bin(dur > .15 & dur < .5,:) = 2;
        bin(dur > .5,:) = 3;
    elseif f == 7
        bin(dur < .2,:) = 1;
        bin(dur > .2 & dur < .5,:) = 2;
        bin(dur > .5,:) = 3;
    elseif f == 10
        bin(dur < .2,:) = 1;
        bin(dur > .2 & dur < .44,:) = 2;
        bin(dur > .44,:) = 3;
    end
elseif strcmp(subject,'521846') && strcmp(session,'Task')
    if f == 8 || f == 9
        bin(dur < .2,:) = 1;
        bin(dur > .2 & dur < .44,:) = 2;
        bin(dur > .47,:) = 3;
    elseif f == 11
        bin(dur < .3,:) = 1;
        bin(dur > .3 & dur < .44,:) = 2;
        bin(dur > .47,:) = 3;
    end
elseif strcmp(subject,'015397') && strcmp(session,'Task')
    if f == 1
        bin(dur < .2,:) = 1;
        bin(dur > .2 & dur < .44,:) = 2;
        bin(dur > .47,:) = 3;
    end
elseif strcmp(subject,'663186') && strcmp(session,'Task')
    if f == 7 || f == 9
        bin(dur < .2,:) = 1;
        bin(dur > .2 & dur < .44,:) = 2;
        bin(dur > .47,:) = 3;
    end
elseif strcmp(subject,'503608') && strcmp(session,'Task')
    if f == 1
        bin(dur <= .16,:) = 1;
        bin(dur > .16 & dur < .47,:) = 2;
        bin(dur >= .47,:) = 3;
    elseif f == 2
        bin(dur <= .2,:) = 1;
        bin(dur > .2 & dur < .37,:) = 2;
        bin(dur >= .37,:) = 3;
    elseif f == 3 || f == 4 || f == 5 || f == 7 || f == 9
        bin(dur <= .2,:) = 1;
        bin(dur > .2 & dur < .41,:) = 2;
        bin(dur >= .41,:) = 3;
    elseif f == 6 || f == 8
        bin(dur <= .23,:) = 1;
        bin(dur > .23 & dur < .48,:) = 2;
        bin(dur >= .48,:) = 3;
    elseif f == 10 || f == 11
        bin(dur <= .24,:) = 1;
        bin(dur > .24 & dur < .48,:) = 2;
        bin(dur >= .48,:) = 3;
    end
elseif strcmp(subject,'147947') && strcmp(session,'Task')
    if f == 8
        bin(dur <= .2,:) = 1;
        bin(dur > .2 & dur < .45,:) = 2;
        bin(dur >= .45,:) = 3;
    end
end

[~,idx] = sort(event.dur);
labels = event.labels(idx);

tlabels = cell(size(triggers,1),1);
for i = 1:size(triggers,1)
    tlabels{i,1} = labels{bin(i)}; 
end
triggers.label = tlabels;

% manual corrections
if strcmp(subject,'088674') && strcmp(session,'Task') && f == 4
    triggers.label{38} = 'Transition'; 
elseif strcmp(subject,'263098') && strcmp(session,'Task') && f == 6
    triggers.label{70} = 'Image';
elseif strcmp(subject,'680913') && strcmp(session,'Task') && f == 2
    triggers(14,:) = []; % for some reason the images didn't play back for trial 4 of block 1 (practice = 0)
elseif strcmp(subject,'383991') && strcmp(session,'Task')
    if f == 1
        triggers.label{45} = 'Choice';
    elseif f == 11
        triggers.label{44} = 'Transition'; 
    end
elseif strcmp(subject,'396430') && strcmp(session,'Task')
    if f == 4
        triggers.label{34} = 'Transition';
    elseif f == 5
        triggers.label{51} = 'Choice';
        triggers.label{59} = 'Choice';
    elseif f == 9
        triggers.label{78} = 'Image';
    elseif f == 10
        triggers.label{25} = 'Choice';
    end
elseif strcmp(subject,'015397') && strcmp(session,'Task')
    if f == 1
        triggers.label{1} = 'Choice';
    elseif f == 3
        triggers.label{4} = 'Image';
        triggers.label{76} = 'Image';
        triggers.label{86} = 'Image'; 
    end
elseif strcmp(subject,'503608') && strcmp(session,'Task')
    if f == 1
        triggers.label{31} = 'Choice'; 
    elseif f == 2
        triggers.label{7} = 'Choice';
        triggers.label{19} = 'Choice';
        triggers.label{51} = 'Choice';
    elseif f == 3
        triggers.label{57} = 'Choice';
        triggers.label{73} = 'Choice';
        triggers.label{91} = 'Choice';
    elseif f == 4
        triggers.label{55} = 'Choice';
    elseif f == 5
        triggers.label{7} = 'Choice';
        triggers.label{13} = 'Choice';
        triggers.label{14} = 'Transition';
        triggers.label{20} = 'Transition';
        triggers.label{25} = 'Choice';
        triggers.label{31} = 'Choice';
        triggers.label{61} = 'Choice';
        triggers.label{63} = 'Choice';
    elseif f == 6
        triggers.label{1} = 'Choice';
        triggers.label{36} = 'Transition';
    elseif f == 7
        triggers.label{7} = 'Choice';
        triggers.label{19} = 'Choice';
        triggers.label{45} = 'Choice';
        triggers.label{55} = 'Choice';
    elseif f == 8
        triggers.label{55} = 'Choice';
    elseif f == 9
        triggers.label{19} = 'Choice';
    elseif f == 11
        triggers.label{1} = 'Choice';
    end
elseif strcmp(subject,'506559') && strcmp(session,'Task')
    if f == 10
        triggers.label{23} = 'Image';
        triggers.label{31} = 'Image';
        triggers.label{58} = 'Image';
        triggers.label{69} = 'Image';
        triggers.label{71} = 'Image';
        triggers.label{74} = 'Image';
    end
end

%% Match to behavioural file

if strcmp(session,'FL')

    bidx = find(event.table.Block == fl.block_idx(f,1));
    B = event.table(bidx,:);

    % see how many trials are in the behavioural log
    bTrials = size(B,1);

    % see how many image/cue presentations are in the MEG file
    mImgIdx = strfindcell(triggers.label,'Image');
    mImgs = size(mImgIdx,1);
    mCueIdx = strfindcell(triggers.label,'Cue');
    mCues = size(mCueIdx,1);

    if bTrials ~= mImgs || bTrials ~= mCues
         error([subject ', ' session ' file ' num2str(f) ': no. of trials does not match'])
    end

    % compare against the behavioural log (RTs)
    for i = 1:mImgs
        bRT = B.RT(i)/1000;
        if i < mImgs
            window = [triggers.onset(mCueIdx(i)) triggers.onset(mCueIdx(i)+1)];
        else
            window = [triggers.onset(mCueIdx(i)) length(trigger)];
        end
        idx = find(responses.onset > window(1) & responses.onset < window(2));
        if length(idx) > 1 % more than one response
            tmp = (responses.onset(idx) - window(1)) / data.fsample;
            idx = idx(find(abs(tmp-bRT) == min(abs(tmp-bRT))));
        end
        mRT = (responses.onset(idx) - window(1)) / data.fsample; % RT recorded by MEG, in seconds
        if abs(bRT-mRT) > .1 % if there's more than a 100ms discrepancy
            error([subject ', ' session ' file ' num2str(f) ' choice ' num2str(i) ' RTs do not match up'])
        end
    end

    % if everything checked out, add labels
    oldtriggers = triggers;
    for i = 1:mImgs
        triggers.label{mImgIdx(i)} = event.table.Image(bidx(i));
        triggers.trial(mImgIdx(i)) = bidx(i);
    end
    for i = 1:mCues
        triggers.label{mCueIdx(i)} = event.table.Acc(bidx(i));
        triggers.trial(mCueIdx(i)) = bidx(i);
    end

elseif strcmp(session,'Task')

    bidx = find(event.table.Block == fl.block_idx(f,1) & event.table.Practice == fl.block_idx(f,2));
    B = event.table(bidx,:);

    if B.Practice == 1
        rtlim = Inf;
    else
        rtlim = 30;
    end

    % see how many choices are in the meg file
    mChoiceIdx = strfindcell(triggers.label,'Choice');
    mChoices = length(mChoiceIdx);

    % compare against the behavioural log
    bChoices = size(B,1);

    if strcmp(subject,'680913') && (f == 2 || f == 3)
        %  several trials (at the beginning of blocks) missing from
        %  behavioural log
    else
        if mChoices ~= bChoices
            error([subject ', ' session ' file ' num2str(f) ': no. of choices does not match'])
        end
    end

    % compare responses (type & speed)   
    bRTs = B.RT;
    bChoices = B.Choice;
    rOnsets = responses.onset;
    if strcmp(subject,'680913') && strcmp(session,'Task')
        if f == 2
            % skip first response because first choice trigger is missing (late recording state)
            bRTs = B.RT(2:end);
            bChoices = B.Choice(2:end);
            rOnsets = responses.onset(2:end);
        end
    end
    for i = 1:mChoices
        window = [triggers.onset(mChoiceIdx(i)) triggers.onset(mChoiceIdx(i)+1)];
        idx = find(rOnsets > window(1) & rOnsets < window(2));
        if strcmp(subject,'680913') && strcmp(session,'Task') && f == 2 && i == 3
            % trial 4 is missing from the behavioural log
        else
            if ~isempty(idx) && bRTs(i) < rtlim
                if length(idx) > 1 % more than one response
                    tmp = [];
                    keepgoing = true;
                    while keepgoing
                        if (rOnsets(idx(1)) - window(1)) / data.fsample < 5
                            idx(1,:) = [];
                        else
                            idx = idx(1,:);
                            keepgoing = false;
                        end
                    end
                end
                mRT = (rOnsets(idx) - window(1)) / data.fsample; % RT recorded by MEG, in seconds
                if ~isnan(bRTs(i))
                    bRTs(i) = mRT;
                else
                    bRT = bRTs(i);
                    if abs(bRT-mRT) > .1 % if there's more than a 100ms discrepancy
                        error([subject ', ' session ' file ' num2str(f) ' choice ' num2str(i) ' RTs do not match up'])
                    end
                    if (rOnsets(idx) == 4 && bChoices(i) ~= 1) || (rOnsets(idx) == 1 && bChoices(i) ~= 2) % response trigger 4 = left (airlock), 1 = right (supply)
                        error([subject ', ' session ' file ' num2str(f) ' choice ' num2str(i) ' type does not match up'])
                    end
                end
            end
        end
    end

    % if everything went through, relabel the triggers
    oldTriggers = triggers;
    triggers.trial = nan(size(triggers,1),1);

    if strcmp(subject,'680913') && strcmp(session,'Task') && f == 2
        % manually enter in the first choice trial (trigger is missing due to late recording start)
        triggers.label{1} = ['Door ' num2str(B.Transition(1))];
        for i = 1:3
            triggers.label{i+1} = event.stimuli.img{bidx(1)}{i};
        end
        triggers.label{5} = 'Outcome';
        triggers.trial(1:5) = B.ExpTrial(1);
    end

    for i = 1:mChoices

        % subset the 'triggers' table to be just this choice trial
        if i < mChoices
            idx = mChoiceIdx(i):mChoiceIdx(i+1)-1;
        else
            idx = mChoiceIdx(i):size(triggers,1);
        end

        % check that the right no. of things happened for this trial type
        pass = true;
        if ~isnan(B.nCombo(i)) % some behavioural trials missing from the log for 680913
            % some behavioural trials missing from the log for 680913
            if bChoices(i) == 1 % risky choice 
                if bRTs(i) < rtlim &&  ~(strcmp(subject,'680913') && strcmp(session,'Task') && f == 2 && i == 3) % made in time
                    rows = {'Choice','Transition','Image','Image','Image','Image'}'; % correct rows
                elseif bRTs(i) >= rtlim || (strcmp(subject,'680913') && strcmp(session,'Task') && f == 2 && i == 3) % made too late
                    rows = {'Choice'}'; % correct rows
                end
                if length(idx) ~= length(rows)
                    if strcmp(subject,'680913') && strcmp(session,'Task') && f == 7 && i == 5
                        % for some reason the image triggers are missing from ExpTrial 107??
                    elseif strcmp(subject,'396430') && strcmp(session,'Task') && f == 6 && i == 18
                        % last trigger got cut off
                    else
                        pass = false;
                    end
                else
                    for r = 1:length(rows)
                        if ~strcmp(triggers.label(idx(r)),rows{r})
                            pass = false;
                        end
                    end
                end
            elseif bChoices(i) == 2 % safe choice
                if bRTs(i) < rtlim % made in time
                    rows = {'Choice','Image'}'; % correct rows
                else % made too late
                    rows = {'Choice'}'; % correct rows
                end
                if length(idx) ~= length(rows)
                    pass = false;
                else
                    for r = 1:length(rows)
                        if ~strcmp(triggers.label(idx(r)),rows{r})
                            pass = false;
                        end
                    end
                end
            else
                warning('');
            end
        end

        if ~pass
            error([subject ', ' session ', file ' num2str(f) ': error with choice ' num2str(i)])
        end

        if isnan(B.nCombo(i))
            if ~isnan(bChoices(i))
                choiceIdx = idx(strfindcell(triggers.label(idx),'Choice'));
                if bChoices(i) == 1
                    choiceLabel = 'Risky choice';
                else choiceLabel = 'Safe choice';
                end
                triggers.label{choiceIdx} = choiceLabel;

                tmp = setdiff(idx,choiceIdx);
                for jj = 1:length(tmp)
                    triggers.label{tmp(jj)} = NaN;
                end
            else
                for jj = 1:length(idx)
                    triggers.label{idx(jj)} = NaN;
                end
            end
        else

            if strcmp(subject,'680913') && strcmp(session,'Task') && f == 7 && i == 5
                triggers.label{11} = NaN;
                triggers.label{12} = NaN;
            else

                % assign the correct choice label (Risky or Safe)
                choiceIdx = idx(strfindcell(triggers.label(idx),'Choice'));
                if bChoices(i) == 1
                    choiceLabel = 'Risky choice';
                else choiceLabel = 'Safe choice';
                end
                triggers.label{choiceIdx} = choiceLabel;

                % assign the correct transition label (Door 1, Door 2, Too slow)
                transIdx = idx(strfindcell(triggers.label(idx),'Transition'));
                if ~isempty(transIdx)
                    if bRTs(i) <= rtlim
                        if strcmp(subject,'680913') && strcmp(session,'Task') && f == 2
                            transLabel = ['Door ' num2str(B.Transition(i+1))];
                        else
                            transLabel = ['Door ' num2str(B.Transition(i))];
                        end
                    else
                        transLabel = 'Too slow';
                    end
                    triggers.label{transIdx} = transLabel;
                end

                % assign the correct image label (state name OR supply room
                imgIdx = idx(strfindcell(triggers.label(idx),'Image'));
                if ~isempty(imgIdx)
                    if bChoices(i) == 2 && length(imgIdx) == 1
                        triggers.label{imgIdx} = 'Supply room';
                    else
                        for j = 1:length(imgIdx)
                            if j < length(imgIdx)
                                if size(event.stimuli.img{bidx(i)},1) == 1
                                    tmp = event.stimuli.img{bidx(i)};
                                    try
                                        if isempty(tmp) || sum(isnan(tmp{1}) == 0)
                                            imgLabel = NaN;
                                        end
                                    catch
                                        imgLabel = NaN;
                                    end
                                else
                                    imgLabel = event.stimuli.img{bidx(i)}{j};
                                end
                            else imgLabel = 'Outcome';
                            end
                            triggers.label{imgIdx(j)} = imgLabel;
                        end
                    end
                end
            end
        end

        % add ID of this trial in the 'event' variable
        if strcmp(subject,'680913') && strcmp(session,'Task') && f == 2
            triggers.trial(idx) = ones(length(idx),1)*B.ExpTrial(i+1);
        else
            triggers.trial(idx) = ones(length(idx),1)*B.ExpTrial(i);
        end

    end
end

fsave = [subject '_' session '_triggers-' num2str(Fs) 'Hz.mat'];

close all 

end
