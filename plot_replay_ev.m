%% Plot EV curves for figures

d = readtable('D:\2020_RiskyReplay\results\ev_replay_differential.csv');

d.Replay_cat = cell(size(d,1),1);
d.Replay_cat(d.Replay==min(d.Replay)) = {'aversive'};
d.Replay_cat(d.Replay==max(d.Replay)) = {'rewarding'};
d.Replay_cat(cellfun('isempty',d.Replay_cat)) = {'ignore'};

mcev = -0.0799505; % what an EV of 1 is after it's been mean-centered

% Plot rewarding & aversive replay
cmap = [14, 255, 127;
    255, 14, 84 ]/255;

figure
ycat = {'diminished','enhanced'};
for i = 1:2
    figure
    if i==1
        tmp=R;
        thiscmap = [225, 255, 239
                    0, 190, 115 ]/255;
    else
        tmp=A;
        thiscmap = [255, 226, 234
                    249, 0, 70 ]/255;
    end
    x = [];
    y = [];
    for c = 1:length(ycat)
        idx = contains(tmp.Replay_cat,ycat{c});
        x = [x, tmp.EV(idx)];
        y = [y, tmp.Choice(idx)];
    end
    
    plot(x(:,1),y(:,1),'color',cmap(i,:),'linewidth',1.4,'linestyle','--'); hold on
    plot(x(:,2),y(:,2),'color',cmap(i,:),'linewidth',1.4); hold on

%     plot([mcev mcev],[0 1],'k:','linewidth',1.5);
    set(gca,'ticklength',[0 0])
end

%% Anxiety & Risk aversion

R = readtable('D:\2020_RiskyReplay\results\ev_replay_PCA1.csv');

R.Replay_cat = cell(size(R,1),1);
R.Replay_cat(R.Replay_differential==min(R.Replay_differential)) = {'aversive'};
R.Replay_cat(R.Replay_differential==max(R.Replay_differential)) = {'rewarding'};
R.Replay_cat(cellfun('isempty',R.Replay_cat)) = {'ignore'};

R.Riskaversion = nan(size(R,1),1);
R.Riskaversion(R.PCA1==min(R.PCA1)) = 1; % high
R.Riskaversion(R.PCA1==max(R.PCA1)) = 0; % low

A = readtable('D:\2020_RiskyReplay\results\ev_replay_PCA2.csv');

A.Replay_cat = cell(size(A,1),1);
A.Replay_cat(A.Replay_differential==min(A.Replay_differential)) = {'aversive'};
A.Replay_cat(A.Replay_differential==max(A.Replay_differential)) = {'rewarding'};
A.Replay_cat(cellfun('isempty',A.Replay_cat)) = {'ignore'};

A.Anxiety = nan(size(A,1),1);
A.Anxiety(A.PCA2==min(A.PCA2)) = 1; % high
A.Anxiety(A.PCA2==max(A.PCA2)) = 0; % low



mcev = -0.0799505; % what an EV of 1 is after it's been mean-centered

% Plot
figure
cc = 0;
for i = 1:2 % anxiety, risk aversion
    if i==1
        tmp = A;
        Y = table2array(tmp(:,contains(tmp.Properties.VariableNames,'Anxiety')));
    else
        tmp = R;
        Y = table2array(tmp(:,contains(tmp.Properties.VariableNames,'Riskaversion')));
    end
    for g = 1:2 % low, high
        cc = cc+1;
        subplot(2,2,cc);
        for c = 1:2 % negative differential replay, positive differential replay
            if c==1
                idx = contains(tmp.Replay_cat,'aversive') & Y==(g-1);
            else
                idx = contains(tmp.Replay_cat,'rewarding') & Y==(g-1);
            end
            x = tmp.EV(idx);
            y = tmp.Choice(idx);
            if c==1
                plot(x,y,'k','linewidth',1.5,'linestyle','--'); hold on
            else
                plot(x,y,'k','linewidth',1.5); hold on
            end
        end
        set(gca,'ticklength',[0 0])
        if g==1
            title('Low');
        else
            title('High');
        end
    end
    if i==1
        sgtitle('Anxiety')
    else
        sgtitle('Risk Aversion')
    end
end