function [output] = ss_plot(d,opts)

% d = trials/subjects x perms x lags
gNames = {'Forward','Backward','Forward-Backward'};

if ~isfield(opts,'x')
    opts.x = 1:size(d,3); % x axis
end
x = opts.x;

if ~isfield(opts,'g')
    opts.g = 3; % 1 = fwd, 2 = bwd, 3 = fwd-bwd
    gNames = {'','','unknown'};
end
g = opts.g;
ttext = gNames{g};

if ~isfield(opts,'showOverall')
    opts.showOverall = true; % whether to show average of main permutation
end
if ~isfield(opts,'showIndividual')
    opts.showIndividual = false; % whether to show each individual traces for main permutation
end
if ~isfield(opts,'showError')
    opts.showError = true;
end

if ~isfield(opts,'nullThresh')
    opts.nullThresh = true;
end
if ~isfield(opts,'nullMoving')
    opts.nullMoving = false;
end
if ~isfield(opts,'subtractNull')
    opts.subtractNull = false;
end

if ~isfield(opts,'nullPerms_overall')
    opts.nullPerms_overall = false; % whether to show the average of each null permutation
end
if ~isfield(opts,'nullPerms_individual')
    opts.nullPerms_individual = false; % whether to show the average of each null permutation
end

if ~isfield(opts,'avCol')
    opts.avCol = [0 0 0]; % colour for average trace
end
avCol = opts.avCol;
if ~isfield(opts,'indCol')
    opts.indCol = colours(size(d,1),'rainbow'); % colour map for individual traces
end
indCol = opts.indCol;

if ~isfield(opts,'linewidth')
    opts.linewidth = 1.2;
end
linewidth = opts.linewidth;

if ~isfield(opts,'twosided')
    opts.twosided = false;
end

if ~isfield(opts,'makeplot')
    opts.makeplot = true;
end

%% Get traces

% Calculate sequenceness mean
y = squeeze(d(:,1,:)); % get first (main) permutation

m = nanmean(y);
sem = squeeze(nanstd(y)/sqrt(size(y,1)));
ul = m + sem;
ll = m - sem;

% Calcualte null
if size(d,2) > 1

    np = squeeze(nanmean(d(:,2:end,:))); % get null permutations, averaged over trials/subjects
    np = unique(np,'rows');

    if g==3 || opts.twosided
        abs_np = unique(round(abs(np),4),'rows');
        npThreshQuant = quantile(max(abs_np,[],2),.975);
        npThreshMax = max(abs(np));
    else
        npThreshQuant = quantile(max(np,[],2),.95);
        npThreshMax = max(np);
    end
else
    npThreshMax = NaN;
    npThreshQuant = nan(size(d,3),1);
end

if opts.subtractNull
    for i = 1:size(y,1)
        npi = squeeze(d(i,2:end,:));
        if g==3 || opts.twosided
            npiThresh = quantile(max(abs(npi),[],2),.975);
            y(i,:) = abs(y(i,:)) - npiThresh;
        else
            npiThresh = quantile(max(npi,[],2),.95);
            y(i,:) = y(i,:) - npiThresh;
        end
    end
%     y = y - npThreshQuant;
    if g == 3 || opts.twosided
        m = abs(m) - npThreshQuant;
    else
        m = m - npThreshQuant;
    end
    ul = m + sem;
    ll = m - sem;
end

output = [];
output.m = m;
output.y = y;
output.ul = ul;
output.ll = ll;
output.npThreshMax = npThreshMax;
output.npThreshQuant = npThreshQuant;

%% Plot

if opts.makeplot
    
    % Show each null permutation (averaged across trials/subjects)
    if opts.nullPerms_overall && size(np,2) > 1
        cmap = colours(size(np,1),'rainbow');
        for perm = 1:size(np,1)
            plot(x,np(perm,:),'color',cmap(perm,:),'linestyle','-'); hold on
        end
    end

    % Plot average sequenceness
    if opts.showOverall
        if opts.showError
            patch([x fliplr(x)],[ul fliplr(ll)],avCol,'facealpha',.2,'edgealpha',0,'handlevisibility','off'); hold on
        end
        plot(x,m,'Color',avCol,'linewidth',linewidth); hold on
    end

    % Plot subject sequenceness
    if opts.showIndividual
        for i = 1:size(y,1)
            plot(x,y(i,:),'color',indCol(i,:),'linewidth',linewidth); hold on
        end
    end

    % Plot subject's maximum null permutation
    if opts.nullPerms_individual
        for i = 1:size(d,1)
            snp = squeeze(d(i,:,:));
            [mn,idx] = max(max(snp,[],2));
            plot(x,snp(idx,:),'color',indCol(i,:),'linewidth',.8*linewidth,'linestyle','--'); hold on
        end
    end

    % Plot null threshold max
    if opts.nullThresh && ~opts.subtractNull
        plot(x([1 end]),[npThreshQuant npThreshQuant],'Color',avCol,'linestyle','--','handlevisibility','off'); hold on
        if g==3 || opts.twosided
            plot(x([1 end]),-[npThreshQuant npThreshQuant],'Color',avCol,'linestyle','--','handlevisibility','off'); hold on
        end
    end

    % Plot moving null threshold
    if opts.nullMoving && ~opts.subtractNull
        plot(x,npThreshMax,'Color',avCol,'linestyle',':','handlevisibility','off'); hold on
        if g==3 || opts.twosided
            plot(x,-npThreshMax,'Color',avCol,'linestyle',':','handlevisibility','off'); hold on
        end
    end

    % Plot zero line
    plot(x([1 end]),[0 0],'k:','handlevisibility','off'); hold on

    % Tidy up axes
    title(gNames{g})
    set(gca,'ticklength',[0 0])
    xlabel('Lags (ms)')
    ylabel('Sequenceness')
end
end