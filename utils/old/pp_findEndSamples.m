function [endSamples, fsave] = pp_findEndSamples(fl,in)

% Manually identify the end of the recording
% (requires FieldTrip has been added to path)

% subject  = character (e.g., '088674')
% fl       = filelist (produced by pp_cfg(subject,session))
% in       = info (produced by pp_cfg(subject,session))

%% Load data

% Get sampling rate
Fs = nan(1,length(fl.data));
for f = 1:length(fl.data)
    hdr = ft_read_header(fullfile(fl.dir_meg,fl.data{f}));
    Fs(f) = hdr.Fs;
end

Fs = unique(Fs);

%% Run manual, graphical check for end of recording

endSamples = [];
for f = 1:length(fl.data)

    % load raw data
    cfg = [];
    cfg.dataset = fullfile(fl.dir_meg,fl.data{f});
    raw = ft_preprocessing(cfg);

    % plot trigger channel
    tC = strfindcell(raw.label,in.triggerChannel);
    trigger = raw.trial{1}(tC,:);
    figure
    plot(trigger)
    title([fl.subject ', file ' num2str(f) ' of ' num2str(length(fl.data))])

    % make interactable
    disp('Click on the end of the recording, then press SPACEBAR')
    set(gcf,'CurrentCharacter',char(1));
    h = datacursormode;
    set(h,'DisplayStyle','datatip','SnapToData','off');
    waitfor(gcf,'CurrentCharacter',char(32));
    tmp = getCursorInfo(h);
    bC = tmp.Position;

    endSamples(f) = bC(1);

end

fsave = [fl.subject '_' fl.session '_endsamples-' num2str(Fs) 'Hz.mat'];

end