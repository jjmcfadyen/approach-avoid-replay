%% Path replay by type & choice

D = readtable('D:\2020_RiskyReplay\results\replay\replay_path_value_choice.csv');

cc = 0;
X = [];
Y = [];
for choice = [1 0]
    for p = 1:2 % rewarding, aversive
        if p==1
            idx = contains(D.Replay_type,'rewarding');
        else
            idx = contains(D.Replay_type,'aversive');
        end
        cc = cc+1;
        Y = [Y, D.Sequenceness(contains(D.Choice,num2str(choice)) & idx)];
        X = [X, ones(size(Y,1),1)*cc];
    end
end
N = size(X,1);

figure
set(gcf,'position',[440 319 365 479])

markersize = 30;

cmap = [0 1 0
        1 0 0];
cmap = repmat(cmap,2,1);

for c = 1:size(X,2)
    [x,y] = beeswarm(Y(:,c),.03,.2);
    X(:,c) = x+c;
end

% (subject lines)
for s = 1:N
    plot(X(s,1:2),Y(s,1:2),'color',[0 0 0 .25]); hold on 
    plot(X(s,3:4),Y(s,3:4),'color',[0 0 0 .25]); hold on 
end

% (subject dots)
for c = 1:size(Y,2)
    scatter(X(:,c),Y(:,c),markersize,'markerfacecolor',cmap(c,:),...
        'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
end

% (boxplot)
width = 0.25;
plot([1+width/2 2-width/2],[median(Y(:,1)) median(Y(:,2))],...
    'k','linewidth',1.6); hold on
plot([3+width/2 4-width/2],[median(Y(:,3)) median(Y(:,4))],...
    'k','linewidth',1.6); hold on
for c = 1:size(Y,2)
    m = median(Y(:,c));
    q = quantile(Y(:,c),[.25 .75]);
    e = quantile(Y(:,c),[.05 .95]);
    patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
        'w','facealpha',.75,'edgecolor','k'); hold on
    plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
end

% (axes)
xlim([0.5 4.5])
set(gca,'box','off')
set(gca,'ticklength',[0 0])
xticks(1:4);
xticklabels({'Approach','Avoid','Approach','Avoid'})
ylabel('Replay')

%% Path temporal modulation by type & choice

D = readtable('D:\2020_RiskyReplay\results\replay\temporalmodulation_choice.csv');

cc = 0;
X = [];
Y = [];
for p = 1:2 % rewarding, aversive
    for choice = [1 0]
        cc = cc+1;
        if p==1
            Y = [Y, D.Modulation_rewarding(D.Choice==choice)];
        else
            Y = [Y, D.Modulation_aversive(D.Choice==choice)];
        end
        X = [X, ones(size(Y,1),1)*cc];
    end
end
N = size(X,1);

figure
set(gcf,'position',[440 319 365 479])

markersize = 30;

cmap = [255, 182, 54
        169, 116, 255]/255;
cmap = repmat(cmap,2,1);

for c = 1:size(X,2)
    [x,y] = beeswarm(Y(:,c),.03,.2);
    X(:,c) = x+c;
end

% (subject lines)
for s = 1:N
    plot(X(s,1:2),Y(s,1:2),'color',[0 0 0 .25]); hold on 
    plot(X(s,3:4),Y(s,3:4),'color',[0 0 0 .25]); hold on 
end

% (subject dots)
for c = 1:size(Y,2)
    scatter(X(:,c),Y(:,c),markersize,'markerfacecolor',cmap(c,:),...
        'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
end

% (boxplot)
width = 0.25;
plot([1+width/2 2-width/2],[median(Y(:,1)) median(Y(:,2))],...
    'k','linewidth',1.6); hold on
plot([3+width/2 4-width/2],[median(Y(:,3)) median(Y(:,4))],...
    'k','linewidth',1.6); hold on
for c = 1:size(Y,2)
    m = median(Y(:,c));
    q = quantile(Y(:,c),[.25 .75]);
    e = quantile(Y(:,c),[.05 .95]);
    patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
        'w','facealpha',.75,'edgecolor','k'); hold on
    plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
end

% (axes)
xlim([0.5 4.5])
set(gca,'box','off')
set(gca,'ticklength',[0 0])
xticks(1:4);
xticklabels({'Approach','Avoid','Approach','Avoid'})
ylabel('Replay')

%% Replay by choice

D = readtable('D:\2020_RiskyReplay\results\replay\choice_replay.csv');
N = length(unique(D.Subject));

Y = [D.Replay(contains(D.Choice,'1')) D.Replay(contains(D.Choice,'0'))];
X = [zeros(size(Y,1),1) ones(size(Y,1),1)];


figure
set(gcf,'position',[440 319 365 479])

markersize = 30;

cmap = [255, 182, 54
        169, 116, 255]/255;

for c = 1:2
    [x,y] = beeswarm(Y(:,c),.03,.2);
    X(:,c) = x+c;
end

% (subject lines)
for s = 1:N
    plot(X(s,:),Y(s,:),'color',[0 0 0 .25]); hold on 
end

% (subject dots)
for c = 1:2
    scatter(X(:,c),Y(:,c),markersize,'markerfacecolor',cmap(c,:),...
        'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
end

% (boxplot) / (EMMs)
width = 0.25;
% plot([1+width/2 2-width/2],[median(Y(:,1)) median(Y(:,2))],...
%     'k','linewidth',1.6); hold on
emm = [-0.00328 0.00231 -0.007809   0.00124  % mean, sem, lower CI, upper CI
       0.00417 0.00242 -0.000586   0.00892];
for c = 1:2
%     m = median(Y(:,c));
%     q = quantile(Y(:,c),[.25 .75]);
    m = emm(c,1);
    q = emm(c,3:4);
    patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
        'w','facealpha',.75,'edgecolor','k'); hold on
    plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
%     plot([c c],[emm(c,1)-emm(c,2) emm(c,1)+emm(c,2)],'k','linewidth',1.2); hold on
%     scatter(c,emm(c,1),markersize*6,'markerfacecolor','k'); hold on
end
plot([1 2],emm(:,1),'k','linewidth',1.6); hold on

% (axes)
xlim([0.5 2.5])
set(gca,'box','off')
set(gca,'ticklength',[0 0])
xticks(1:2);
xticklabels({'Approach','Avoid'})
ylabel('Differential Replay')

%% Replay by choice AND personality (median split)

D = readtable('D:\2020_RiskyReplay\results\replay\choice_replay_personality.csv');

for p = 1:2 % anxiety, risk aversion
    
    if p==1
        idx = contains(D.HighAnx,'TRUE');
    else
        idx = contains(D.HighRiskAverse,'TRUE');
    end
    
    cc = 0;
    X = [];
    Y = [];
    for choice = [1 0]
        for g = [0 1]
            cc = cc+1;
            Y = [Y, D.Replay(contains(D.Choice,num2str(choice)) & idx==g)];
            X = [X, ones(size(Y,1),1)*cc];
        end
    end
    N = size(X,1);

    figure
    set(gcf,'position',[440 319 365 479])

    markersize = 30;

    cmap = [255, 182, 54
            169, 116, 255]/255;
    cmap = repmat(cmap,2,1);

    for c = 1:size(X,2)
        [x,y] = beeswarm(Y(:,c),.03,.2);
        X(:,c) = x+c;
    end

    % (subject lines)
    for s = 1:N
        plot(X(s,1:2),Y(s,1:2),'color',[0 0 0 .25]); hold on 
        plot(X(s,3:4),Y(s,3:4),'color',[0 0 0 .25]); hold on 
    end

    % (subject dots)
    for c = 1:size(Y,2)
        scatter(X(:,c),Y(:,c),markersize,'markerfacecolor',cmap(c,:),...
            'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
    end

    % (boxplot)
    width = 0.25;
    plot([1+width/2 2-width/2],[median(Y(:,1)) median(Y(:,2))],...
        'k','linewidth',1.6); hold on
    plot([3+width/2 4-width/2],[median(Y(:,3)) median(Y(:,4))],...
        'k','linewidth',1.6); hold on
    for c = 1:size(Y,2)
        m = median(Y(:,c));
        q = quantile(Y(:,c),[.25 .75]);
        e = quantile(Y(:,c),[.05 .95]);
        patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
            'w','facealpha',.75,'edgecolor','k'); hold on
        plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
    end

    % (axes)
    xlim([0.5 4.5])
    set(gca,'box','off')
    set(gca,'ticklength',[0 0])
    xticks(1:4);
    xticklabels({'Approach','Avoid','Approach','Avoid'})
    ylabel('Differential Replay')
    if p==1
        sgtitle('Anxiety')
    elseif p==2
        sgtitle('Risk Aversion')
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap = [255, 0, 93
    0, 242, 154]/255;

errortype = 'sem'; % 'ci' or 'sem'

figure
for i = 1:4

    if i < 3
        personality = 'anx';
    else
        personality = 'risk';
    end
    if i==1 || i==3
        plevel = 'low';
    else
        plevel = 'high';
    end
    
    D = readtable(['D:\2020_RiskyReplay\results\replay\ev_replay_' personality '-' plevel '.csv']);
    
    if strcmp(errortype,'sem')
        sem = (D.Aversive_ymax-D.Aversive_ymin)/3.92;
        D.Aversive_ymax = D.Aversive_Choice+sem;
        D.Aversive_ymin = D.Aversive_Choice-sem;
        
        sem = (D.Rewarding_ymax-D.Rewarding_ymin)/3.92;
        D.Rewarding_ymax = D.Rewarding_Choice+sem;
        D.Rewarding_ymin = D.Rewarding_Choice-sem;
    end
    
    x = D.EV;
    
    subplot(2,2,i)
    
%     % plot neither bias (black dotted line)
%     plot(x,D.Neither_Choice,'k:','linewidth',1.4); hold on
    
    % plot ribbons
    patch([x' fliplr(x')],[D.Aversive_ymin' fliplr(D.Aversive_ymax')],cmap(1,:),'facealpha',.2,'edgealpha',0); hold on
    patch([x' fliplr(x')],[D.Rewarding_ymin' fliplr(D.Rewarding_ymax')],cmap(2,:),'facealpha',.2,'edgealpha',0); hold on
    
    % plot lines
    plot(x,D.Aversive_Choice,'color',cmap(1,:),'linewidth',1.3); hold on
    plot(x,D.Rewarding_Choice,'color',cmap(2,:),'linewidth',1.3); hold on
    
    set(gca,'ticklength',[0 0])
    box('off')
    
end
