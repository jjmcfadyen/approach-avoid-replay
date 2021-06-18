function [h] = emmip(emm, groupformula, varargin)
% Plotting estimated marginal means (only for categorical vars)
% Copyright 2019, J. Hartman, Cornell University
% written using Matlab version R2016a 
% (errorbar,linkaxes behavior or output changes in later versions)
% 
% INPUTS:
% emm is output from emmeans.m
% groupformula contains variable names for grouping to compare (VAR1), to 
% show dependence (VAR2) and to plot together (VAR3). The input structure 
% mirrors emmip in R, with e.g. VAR1 ~ VAR2 | VAR3. Only VAR2 is required.
% varargin is just for future error option ('SE' or 'CI') ['CI' not yet included]
% 
% example:
% emmT = emmeans(mdl);
% h = emmip(emmT, '~ Smoker');
[var2compare, var2indep, var2group, err] = parseInputs(emm,groupformula,varargin);
% only var2indep must be defined
ind_names = unique(emm.table{:,var2indep},'stable');
% design figure based on var2group
if isempty(var2group) % no grouping
    plot_names{1} = '';
    data_table{1} = emm.table;
    ncol = 1;
    nrow = 1;
else % grouping by 1+ variables
    for gg=length(var2group):-1:1
        [gid(:,gg),plot_names{gg}] = grp2idx(emm.table.(var2group{gg}));
    end
    % design figure subplots
    [npn,ia] = sort(cellfun('length', plot_names),'ascend');
    plot_names = plot_names(ia);
    var2group = var2group(ia);
    gid = gid(:,ia);
    if length(var2group)>1 % 2 grouping variables
        nrow = npn(1);
        ncol = npn(2);
    else % 1 grouping variable
        ncol = min(npn,5); % at most 5 columns of subplots
        nrow = ceil(npn/ncol);
    end
    
    igroup = sortrows(unique(gid,'rows','stable'),1:size(gid,2)-1);
    for ig=size(igroup,2):-1:1
        ugroup(:,ig) = plot_names{ig}(igroup(:,ig));
    end
    % there's a warning for ismember('rows') with cell array input, but it
    % still works so i silence it
    warning('off','MATLAB:ISMEMBER:RowsFlagIgnored');
    for ug=size(ugroup,1):-1:1
        data_table{ug} = emm.table( all(ismember(emm.table{:,var2group},ugroup(ug,:),'rows'),2), : );
    end
    warning('on','MATLAB:ISMEMBER:RowsFlagIgnored');
end
% plot data based on var2compare
figure,
nplots = length(data_table);
for np=nplots:-1:1
    if nplots>1
        ax(np) = subplot(nrow,ncol,np);
    else
        ax(np) = gca;
    end
    if isempty(var2compare) % draw single line per plot
        data2plot = data_table{np}.Estimated_Marginal_Mean;
        err2plot = data_table{np}.(err);
        allpnames = cat(1,plot_names{:});
        if nplots>1 && length(var2group)==2
            [ir,ic] = ind2sub([nrow ncol],np); % only works since length(var2group)<=2
            datanames = {strjoin(allpnames([ir,nrow+ic])',{'_'})};
        elseif nplots>1 && length(var2group)==1
            datanames = allpnames(np);
        else
            datanames = plot_names;
        end
    else % draw multiple lines per plot
        ucompare = unique(emm.table(:,var2compare),'stable');
        for uc=height(ucompare):-1:1
            table2plot = data_table{np}( ismember(data_table{np}(:,var2compare),ucompare(uc,:)), {'Estimated_Marginal_Mean',err} );
            data2plot(:,uc) = table2plot.Estimated_Marginal_Mean;
            err2plot(:,uc) = table2plot.(err);
            datanames{uc} = strjoin( ucompare{uc,:},{'\_'} );
        end
    end
    h{np} = errorbar( data2plot, err2plot );
    % modify aesthetics
    aesmod(h{np},datanames); % color, line style, markers and display names
    % modify limits
    ax(np).XLim = [0 length(ind_names)+1];
    ax(np).YLim = [min(min(data2plot - err2plot - 0.5)),...
                   max(max(data2plot + err2plot + 0.5))];
end
%% beautify plot
% link axes 
if ~isempty(regexp(version,'R2016a','once')) % in 2016a linkaxes scales all to 1st axis rather than max/min
    axlims = reshape([ax(:).YLim]',2,[])';
    for al=1:nplots
        ax(al).YLim = [min(axlims(:,1)) max(axlims(:,2))];
    end
end
linkaxes(ax,'xy');
% add titles and axis labels
if isempty(var2group)
    ax(1).Title.String = plot_names{1};
    xlab = reshape([var2indep' repmat({'_'},size(var2indep,2),1)]',1,[]);
    ax(1).XLabel.String = [xlab{1:end-1}];
    ax(1).YLabel.String = 'Predicted values';
    indepplots = 1;
    rightplots = 1;
elseif length(var2group)==1
    for tp=nplots:-1:1
        ax(tp).Title.String = plot_names{1}{tp};
    end
    indepplots = nplots:-1:1;
    rightplots = ncol:ncol:nplots;
else
    for tp=ncol:-1:1
        ax(tp).Title.String = plot_names{2}{tp};
    end
    leftplots = 1:ncol:nplots;
    for lp=length(leftplots):-1:1
        ax(leftplots(lp)).YLabel.String = plot_names{1}{lp};
    end
    indepplots = nplots:-1:nplots-ncol+1;
    rightplots = ncol:ncol:nplots;
end
% only bottom plots have Xticks and Xtick labels
for ip=1:length(indepplots)
    ax(indepplots(ip)).XTickLabelRotation = 60;
    rotationoffset = 0; %0.25;
    ax(indepplots(ip)).XTick = ax(1).XLim(1)-rotationoffset:ax(1).XLim(2)-rotationoffset;
    ax(indepplots(ip)).XTickLabel = {'',ind_names{:},''};
end
otherxplots = setdiff(1:nplots,indepplots);
for oxp=1:length(otherxplots)
    ax(otherxplots(oxp)).XTick = [ax(1).XLim(1) ax(1).XLim(2)];
    ax(otherxplots(oxp)).XTickLabel = {'',''};
end
% only right plots have Yticks and values, and if multiple lines a legend
for rp=1:length(rightplots)
    if ~isempty(var2group)
        ax(rightplots(rp)).YAxisLocation = 'right';
    end
    if ~isempty(var2compare)
        legend(ax(ncol),'show');
    end
end
otheryplots = setdiff(1:nplots,rightplots);
for oyp=1:length(otheryplots)
    ax(otheryplots(oyp)).YTick = [ax(1).YLim(1) ax(1).YLim(2)];
    ax(otheryplots(oyp)).YTickLabel = {''};
end
end
%% helper functions
function [var2compare,var2indep,var2group,err] = parseInputs(emm,grpformula,erropt)
% clean inputs from the formula
grpformula = strrep(grpformula,' ','');
% error check first
allgrps = regexp(grpformula,'[*~|]','split');
allgrps(strcmp(allgrps,'')) = []; % if no variable on one side of separator
if any(~ismember(allgrps,emm.table.Properties.VariableNames))
    error('Oops. Inputs to emmip are not all valid variables or an incorrect separator was used.');
end
if any(~ismember(emm.table.Properties.VariableNames(1:end-3),allgrps))
    error('All valid variables must be included in the formula for a given EMM input.');
end
if length(allgrps) ~= length(unique(allgrps))
    error('Variables cannot be included more than once in the formula.');
end
if length(strfind(grpformula,'~'))>1 || length(strfind(grpformula,'|'))>1
    error('Error in formula. Too many separators.');
end
if ~ismember(erropt,{'SE', 'CI'})
    error('Can only plot ''SE'' or ''CI''.');
end
% defaults
err = 'SE';
% unpack formula more thoroughly
if isempty(regexp(grpformula,'~','once')) % no '~' in formula
    var2compare = [];
    if isempty(regexp(grpformula,'\|','once')) % no '~' nor '|' 
        var2group = [];
        var2indep = regexp(grpformula,'*','split');
    else % no '~' but there is a '|' 
        conds = regexp(grpformula,'\|','split');
        var2group = regexp(conds{2},'*','split');
        var2indep = regexp(conds{1},'*','split');
    end
else % '~' in formula
    squig = regexp(grpformula,'~','split');
    if strcmp(squig{1},'') % '~' but nothing on left side
        var2compare = [];
    else % '~' with variables on both sides
        var2compare = regexp(squig{1},'*','split');
    end
    if isempty(regexp(squig{2},'\|','once')) % '~' but no '|'
        var2group = [];
        var2indep = regexp(squig{2},'*','split');
	else % '~' and there is a '|'
        conds = regexp(squig{2},'\|','split');
        var2group = regexp(conds{2},'*','split');
        var2indep = regexp(conds{1},'*','split');
    end
end
% can't really plot if too many grouping variables
if length(var2group)>2
    error('Cannot make 2D plot with > 2 grouping variables.');
end
% for now, can't plot if too many comparison variables
% (could use different markers or linestyles in future)
if length(var2compare)>2
    error('Cannot have more than 1 comparison variable. Could be updated in future.');
end
% error option
if ~isempty(erropt) && strcmp(erropt,'CI')
    err = emm.table.Properties.VariableNames{end};
end
end
function [] = aesmod(h,dispnames)
nlines = length(h);
lineflag = false;
if nlines==1
    colors = [0 0 0];
    lineflag = true;
elseif nlines<4
    colors = [1 0 0; 0 1 0; 0 0 1];
    lineflag = true;
else
    colors = parula(nlines);
end
for ii=1:nlines
	h(ii).Color = colors(ii,:);
    h(ii).Marker = '.';
    h(ii).MarkerSize = 10;
    if lineflag
        h(ii).LineStyle = ':';
    end
    h(ii).DisplayName = dispnames{ii};
end
    
end
