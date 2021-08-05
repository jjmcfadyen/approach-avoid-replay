results = parse_behav(subjects{22},'D:\2020_RiskyReplay\data\behav');

figure
set(gcf,'position',[680 424 1097 554])
markersize = 25;
markeralpha = 0.8;
linewidth = 0.5;

% %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,1)

plot([1 size(results,1)],[0 0],'k:','linewidth',1.5); hold on

catchtrials = find((results.nV_1>0 & results.nV_2>0) | (results.nV_1<0 & results.nV_2<0));
for i = 1:length(catchtrials)
    plot(repmat(catchtrials(i),2,1),[-12 12],'color',[255, 214, 32]/255,'linewidth',2); hold on 
end

plot(1:size(results,1),results.nV_1,'k'); hold on
scatter(find(results.nV_1>0),results.nV_1(results.nV_1>0),markersize,'markerfacealpha',markeralpha,'markerfacecolor',[0, 240, 160]/255,'markeredgecolor','k','linewidth',markeredge); hold on
scatter(find(results.nV_1<0),results.nV_1(results.nV_1<0),markersize,'markerfacealpha',markeralpha,'markerfacecolor',[255, 0, 93 ]/255,'markeredgecolor','k','linewidth',markeredge); hold on

plot(1:size(results,1),results.nV_2,'k'); hold on
scatter(find(results.nV_2>0),results.nV_2(results.nV_2>0),markersize,'markerfacealpha',markeralpha,'marker','^','markerfacecolor',[0, 240, 160]/255,'markeredgecolor','k','linewidth',markeredge); hold on
scatter(find(results.nV_2<0),results.nV_2(results.nV_2<0),markersize,'markerfacealpha',markeralpha,'marker','^','markerfacecolor',[255, 0, 93 ]/255,'markeredgecolor','k','linewidth',markeredge); hold on

set(gca,'ticklength',[0 0])
xlim([1 size(results,1)])
set(gca,'box','off')
set(gca,'xtick',[])
set(gca,'ytick',[])

% %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,2)

cmap = colours(5,'viridis');
probs = [.1 .3 .5 .7 .9];

plot([1 size(results,1)],[0.5 0.5],'k:','linewidth',1.5); hold on
plot(1:size(results,1),results.P,'k'); hold on

for i = 1:length(probs)
    scatter(find(results.P==probs(i)),results.P(results.P==probs(i)),markersize,'markerfacealpha',markeralpha,'markerfacecolor',cmap(i,:),'markerfacealpha',1,'markeredgecolor','k','linewidth',markeredge); hold on
end

set(gca,'ticklength',[0 0])
xlim([1 size(results,1)])
set(gca,'box','off')
set(gca,'xtick',[])
set(gca,'ytick',[])

% %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,3)

negCombos = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]; % negator combos
negs = zeros(size(results,1),2);
for i = 1:size(results,1)
    negs(i,:) = negCombos(results.nCombo(i),:); 
end

imagesc(negs');
colormap(colours(6,'viridis'))
set(gca,'ticklength',[0 0])
set(gca,'box','off')
set(gca,'xtick',[])
set(gca,'ytick',[])

% %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,4)

plot([1 size(results,1)],[1 1],'k:'); hold on

plot(1:size(results,1),results.EV,'k'); hold on

scatter(find(results.EV>1),results.EV(results.EV>1),markersize,'markerfacealpha',markeralpha,'markerfacecolor',[0, 240, 160]/255,'markeredgecolor','k','linewidth',markeredge); hold on
scatter(find(results.EV<1),results.EV(results.EV<1),markersize,'markerfacealpha',markeralpha,'markerfacecolor',[255, 0, 93 ]/255,'markeredgecolor','k','linewidth',markeredge); hold on
scatter(find(results.EV==1),results.EV(results.EV==1),markersize,'markerfacealpha',markeralpha,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k','linewidth',markeredge); hold on

set(gca,'ticklength',[0 0])
xlim([1 size(results,1)])
set(gca,'box','off')
set(gca,'ytick',[])
