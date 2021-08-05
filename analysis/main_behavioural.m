%% Behaviour

clear all
clc

%% Directories & parameters

addpath('utils');

dir_data = 'D:\2020_RiskyReplay\data\behav';
dir_raw = 'D:\2020_RiskyReplay\data\meg\raw';

parameters = get_parameters(dir_raw);

subjects = unique(parameters.schar);
N = length(subjects);

%% Plot accuracy & RT & EV

acc = nan(N,3);
rt = nan(N,3);
ev = nan(N,2);
tc = nan(N,2);
optimal = nan(N,2);
for s = 1:N
   
    % Read in data
    d = parse_behav(subjects{s},dir_data);
    
    % Exclude missed trials
    d = d(d.RT<30,:);
    
    % Exclude forced-choice trials
    d = d(d.Forced==0,:);
    
    % Log results
    for c = 1:2
        acc(s,c) = mean(d.Acc(d.Choice==c)); 
        rt(s,c) = mean(d.RT(d.Choice==c));
        ev(s,c) = mean(d.EV(d.Choice==c));
        tc(s,c) = mean(d.Choice==c);
    end
    acc(s,3) = mean(d.Acc);
    rt(s,3) = mean(d.RT);
    
    optimal(s,1) = mean(d.EV >= 1);
    optimal(s,2) = mean(d.EV <= 1);
    
end

%% Statistics

excludedSubjects = {'263098','680913'};
includedidx = ~ismember(subjects,excludedSubjects);

disp('-----------------------------------------')
disp(['Pcnt approach = ' num2str(round(mean(tc(includedidx,1)),4))])
disp(['Pcnt avoid = ' num2str(round(mean(tc(includedidx,2)),4))])
disp('-----------------------------------------')

disp('-----------------------------------------')
disp(['Acc overall = ' num2str(round(mean(acc(includedidx,3)),4)) ', sd = ' num2str(round(std(acc(includedidx,3)),4))])
disp(['Acc approach = ' num2str(round(mean(acc(includedidx,1)),4)) ', sd = ' num2str(round(std(acc(includedidx,1)),4))])
disp(['Acc avoid = ' num2str(round(mean(acc(includedidx,2)),4)) ', sd = ' num2str(round(std(acc(includedidx,2)),4))])
disp('-----------------------------------------')

disp('-----------------------------------------')
disp(['RT approach = ' num2str(round(mean(rt(includedidx,1)),3)) ', sd = ' num2str(round(std(rt(includedidx,1)),3))])
disp(['RT avoid = ' num2str(round(mean(rt(includedidx,2)),3)) ', sd = ' num2str(round(std(rt(includedidx,2)),3))])
disp('-----------------------------------------')

disp('-----------------------------------------')
disp(['EV approach = ' num2str(round(mean(ev(includedidx,1)),3)) ', sd = ' num2str(round(std(ev(includedidx,1)),3))])
disp(['EV avoid = ' num2str(round(mean(ev(includedidx,2)),3)) ', sd = ' num2str(round(std(ev(includedidx,2)),3))])
disp('-----------------------------------------')

% Compare choices on accuracy and RT
[~,p1,ci1,stats1] = ttest(acc(includedidx,1),acc(includedidx,2),'tail','both');
[~,p2,ci2,stats2] = ttest(rt(includedidx,1),rt(includedidx,2),'tail','both');
[~,p3,ci3,stats3] = ttest(ev(includedidx,1),ev(includedidx,2),'tail','both');

writetable(array2table([acc rt],'variablenames',{'acc_approach','acc_avoid','rt_approach','rt_avoid'}),...
    fullfile('D:\2020_RiskyReplay\results','accrt.csv'));

T = [];
for c = 1:2
    thist = array2table([[1:N]' acc(:,c) rt(:,c) ones(N,1)*c],...
        'variablenames',{'Subject','Acc','RT','Choice'});
    T = [T; thist];
end
T.Subject = categorical(T.Subject,unique(T.Subject),subjects);
T.Choice = T.Choice+10;
T.Choice(T.Choice==11) = 1;
T.Choice(T.Choice==12) = 0;

glme = fitglme(T,'Acc~Choice*RT+(1|Subject)')

%% Plot

figure
set(gcf,'position',[440 471 769 327])

markersize = 20;

for p = 1:3
    
    subplot(1,3,p)

    cmap = [255, 182, 54
            169, 116, 255]/255;

    X = nan(N,2);
    Y = nan(N,2);
    for c = 1:2
        if p==1
            [x,y] = beeswarm(acc(:,c),.05,.2);
        elseif p==2
            [x,y] = beeswarm(ev(:,c),.2,.2);
        elseif p==3
            [x,y] = beeswarm(rt(:,c),.2,.2);
        end
        X(:,c) = x+c;
        Y(:,c) = y;
    end

    % (subject lines)
    for s = 1:N
        plot(X(s,:),Y(s,:),'color',[0 0 0 .25]); hold on 
    end

    % (subject dots)
    for c = 1:2
        scatter(X(includedidx,c),Y(includedidx,c),markersize,'markerfacecolor',cmap(c,:),...
            'markerfacealpha',.75,'markeredgealpha',1,'markeredgecolor','k'); hold on
        if any(~includedidx)
            scatter(X(~includedidx,c),Y(~includedidx,c),markersize,'marker','x',...
            'markeredgecolor','k','markeredgealpha',.5); hold on
        end
    end

    % (boxplot)
    width = 0.25;
    plot([1+width/2 2-width/2],[median(Y(includedidx,1)) median(Y(includedidx,2))],...
        'k','linewidth',1.6); hold on
    for c = 1:2
        m = median(Y(includedidx,c));
        q = quantile(Y(includedidx,c),[.25 .75]);
        e = quantile(Y(includedidx,c),[.05 .95]);
        patch([c-width/2 c+width/2 c+width/2 c-width/2],[q(2) q(2) q(1) q(1)],...
            'w','facealpha',.75,'edgecolor','k'); hold on
        plot([c-width/2 c+width/2],[m m],'k','linewidth',1.4); hold on
    %     plot([c c],[q(1) e(1)],'k','linewidth',1.2); hold on
    %     plot([c c],[q(2) e(2)],'k','linewidth',1.2); hold on
    end

    % (axes)
    xlim([0 3])
    set(gca,'box','off')
    set(gca,'ticklength',[0 0])
    xticks(1:2);
    xticklabels({'Approach','Avoid'})
    if p==1
        ylabel('Accuracy')
        ylim([0.4 1])
    elseif p==2
        ylabel('EV')
        plot([0 3],[1 1],'k--','linewidth',1.2); hold on
    elseif p==3
        ylabel('Decision Time')
        ylim([5 17])
    end
end
