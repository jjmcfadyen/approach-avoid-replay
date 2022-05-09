function [onsets] = get_replayOnsets(params,classifier,lagrange,threshType,behav)
% Identify replay onsets
% lagrange = in samples

bcwindow = [-.1 0]; % baseline correct planning window

%% Parameters

dir_meg = 'D:\2020_RiskyReplay\data\meg';

subject = params.schar{1};

nStates = 6;

U = generate_nullperms([]); 
nPerms = size(U,1);

%% Get onsets per run

% Create R vector
R = cell(1,length(lagrange));
R_rewarding = cell(1,length(lagrange));
R_aversive = cell(1,length(lagrange));

% Cycle through runs
onsets = [];
for run = 1:size(params,1)

    disp(['Getting onsets for ' subject ', run ' num2str(run) ' of ' num2str(size(params,1)) '...'])

    % Get task data (100 Hz) for this run
    D = spm_eeg_load(fullfile(dir_meg,['5_ICA_ds-100Hz'],subject,['ICA_ds-100Hz_' subject '_task_r' num2str(params.block(run)) '.mat'])); % loads 'epoched' variable
    
    % Get trial segments
    events = D.events;
    eventtypes = extractfield(events,'type');
    eventtimes = extractfield(events,'time');
    eventvals = extractfield(events,'value');
    if iscell(eventvals)
        tmp = eventvals;
        eventvals = nan(length(eventvals),1);
        for i = 1:length(eventvals)
            if isnumeric(tmp{i})
                eventvals(i) = tmp{i};
            else
                eventvals(i) = NaN;
            end
        end
    end
    trlidx = find(contains(eventtypes,'decision'));
    trialtimes = eventtimes(trlidx)';
    trialtimes(:,2) = NaN;
    trialnums = nan(size(trialtimes,1),1);
    for trl = 1:size(trialtimes,1)
        if length(eventtimes) >= trlidx(trl)+1
            trialnums(trl,1) = eventvals(trlidx(trl));
            if strcmp(eventtypes(trlidx(trl)+1),'outcome') || strcmp(eventtypes(trlidx(trl)+1),'transition')
                trialtimes(trl,2) = eventtimes(trlidx(trl)+1);
            end
        end
    end

    % ignore forced-choice trials & trials with no clear offset
    if params.block(run)==0
        ridx = 1:6;
    else
        ridx = 1:4;
    end
    ridx = unique([ridx find(any(isnan(trialtimes')))]);

    trialtimes(ridx,:) = [];
    trialnums(ridx,:) = [];

    nTrls = size(trialtimes,1);

    % Loop through trials
    for trl = 1:nTrls
    
        if trialtimes(trl,2)-trialtimes(trl,1) >= 5 % trial must be at least 5 seconds long
    
            % Extract data plus baseline
            sampidx = [D.indsample(trialtimes(trl,1)+bcwindow(1)):D.indsample(trialtimes(trl,2))];
            chanidx = setdiff(1:size(D,1),[find(~contains(D.chantype,'MEGGRAD'))]);

            X = squeeze(D(chanidx,sampidx,:))'; % channels x samples

            % Baseline-correct
            bc = mean(X(1:round(abs(bcwindow(1))*D.fsample),:));
            X = X - bc;

            X = X(round(abs(bcwindow(1))*D.fsample)+1:end,:);
            sampidx = sampidx(round(abs(bcwindow(1))*D.fsample)+1:end);

            % Scale data
            X = X ./ prctile(abs(X(:)),95);
        
            % Get reactivation matrix
            X = X*classifier.betas';
        
            % Time-shift per lag
            thisR = nan(nPerms,size(X,1),length(lagrange),2);
            for perm = 1:nPerms
                for path = 1:2
                    for lag = 1:length(lagrange)
                
                        Xdelta = [X(lagrange(lag):end,:); zeros(lagrange(lag)-1,nStates)];
                        XdeltaP = [];
                        if path==1
                            XdeltaP = [X(:,U(perm,1)).*Xdelta(:,U(perm,2)) X(:,U(perm,2)).*Xdelta(:,U(perm,3))];
                        else
                            XdeltaP = [X(:,U(perm,4)).*Xdelta(:,U(perm,5)) X(:,U(perm,5)).*Xdelta(:,U(perm,6))];
                        end
                        thisR(perm,:,lag,path) = sum(XdeltaP,2);
                
                    end
                end
            end
        
            % Match this epoch's time period to the full run (ICA file)
            run_samples = sampidx;
            run_times = sampidx * (1/D.fsample);

            % Get path reactivation strength
            reactivation = squeeze(max(thisR,[],4));

            if params.block(run)==0
                thisblock = 1;
                thispractice = 1;
            else
                thisblock = params.block(run);
                thispractice = 0;
            end
            thistrial = sprintf('%04d',trialnums(trl));
            thistrial = str2double(thistrial(end-1:end));
            behavidx = behav.Practice==thispractice & behav.Block==thisblock & behav.Trial==thistrial;
            if sum(behavidx)~=1
                error('Wrong behavioural index')
            end

            pathtypes = [behav.nV_1(behavidx) behav.nV_2(behavidx)];

            if (pathtypes(1)>0 && pathtypes(2)<0) || (pathtypes(1)<0 && pathtypes(2)>0)
                
                rewarding = squeeze(thisR(:,:,:,find(pathtypes>0)));
                aversive = squeeze(thisR(:,:,:,find(pathtypes<0)));

                % Add to variables
                for lag = 1:length(lagrange)
                    R{lag} = [R{lag}; squeeze(reactivation(:,:,lag))'];
                    R_rewarding{lag} = [R_rewarding{lag}; squeeze(rewarding(:,:,lag))'];
                    R_aversive{lag} = [R_aversive{lag}; squeeze(aversive(:,:,lag))'];
                end
                onsets = [onsets;
                    [ones(length(run_samples),1)*params.block(run) run_samples',run_times']
                    ];
            end
        end
    end
end

O = array2table(onsets,'variablenames',{'Run','Sample','Time'});

onsets = [];
for lag = 1:length(lagrange)
    thisO = O;
    thisO.Reactivation = R{lag}(:,1);
    thisO.Path = R_rewarding{lag}(:,1) > R_aversive{lag}(:,1);
    thisO.Lag = ones(size(thisO,1),1)*lagrange(lag)*10;
    onsets = [onsets; thisO];
end

[~,sortidx] = sortrows([onsets.Run onsets.Sample onsets.Lag]);
onsets = onsets(sortidx,:);
    
% concatenate all the lags
Rnull = [];
for lag = 1:length(lagrange)
    Rnull = [Rnull; R{lag}(:,2:end)];
end

switch threshType
    case 'nullperm'
        Rthresh = quantile(max(Rnull,[],2),.95); % get null threshold across all lags (should make best lag win)
    case 'toppercent'
        Rthresh = quantile(onsets.Reactivation,.95);
end
onsets.Onset = onsets.Reactivation > Rthresh;

onsets = onsets(onsets.Onset==true,:);

% remove any onsets where there are onsets in the previous 100 ms
window = [-100 0]; % in ms
d = [NaN; diff(onsets.Time)];
onsets.Onset(d<abs(window(1)/1000)) = false;

onsets = onsets(onsets.Onset==true,:);

end