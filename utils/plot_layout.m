function plot_layout(data,labels)

% D is a fieldtrip object

%% Get layout

cfg = [];
cfg.layout = 'CTF275_helmet.lay';
layout = ft_prepare_layout(cfg);

L = layout.label;

%% Get data

Z = nan(1,length(L));
for chan = 1:length(L)
    idx = strfindcell(labels,L{chan},'find');
    if ~isempty(idx)
        Z(chan) = nanmean(data(strfindcell(labels,L{chan}),:));
    else
        Z(chan) = NaN; 
    end
end

idx = ~isnan(Z);

Z = Z(idx)';
X = layout.pos(idx,1);
Y = layout.pos(idx,2);
L = L(idx);

Z = (Z - min(Z)) / ( max(Z) - min(Z) ); % normalise between 0 and 1
Z = Z + .01;

%% Plot

figure

cmap = colours(101,'viridis');

% plot points
for i = 1:length(Z)
    scatter(X(i), Y(i), 'markeredgealpha',0, 'markerfacecolor', cmap(ceil(Z(i)*100),:)); hold on
end

% plot outline
for i=1:length(layout.outline)
    if ~isempty(layout.outline{i})
        X = layout.outline{i}(:,1);
        Y = layout.outline{i}(:,2);
        plot(X, Y); hold on
    end
end

% clean up
set(gca,'TickLength',[0 0])

end
