
% Run assess_sequenceness.m up until line 144 (specifiying the example subject in the for 
% loop starting on line 62). Go as far as line 144 for that subject, then execute this code.

%% Get sequenceness

this_sq = [];
for trl = 1:nTrls
    
    pred = cc_predict(test.D{trl},classifier,1,1);
    tmp =  ss_build(pred,U,in);
    this_sq(trl,:,:,:,:) = squeeze(tmp(:,:,:,:,:));
    
end

opts = [];
opts.g = 1;
opts.x = lags;
% opts.showIndividual = true;
% opts.subtractNull = true;
figure
for tr = 1:4
    subplot(2,2,tr)    
    ss_plot(squeeze(this_sq(:,:,opts.g,tr,:)),opts); hold on
end

%% Get best trials for both paths

plagrange = [30 80];

trial_evidence = [];
for trl = 1:size(this_sq,1)
   
    pidx = lags >= plagrange(1) & lags <= plagrange(2);
    
    y = squeeze(mean(this_sq,3));
    
    np = squeeze(y(trl,2:end,pidx));
    npthresh = quantile(abs(max(np,[],2)),.95);
    
    athresh = squeeze(y(trl,1,pidx))' - npthresh;
    trial_evidence(trl,:) = max(athresh);
    
end

figure
plot(trial_evidence); hold on

[~,best_trial] = max(trial_evidence');

%% Plot heatmap for specific trial

L = 50;
smoothing = 3;

all_evidence = zeros(length(test.x),2);
for trl = 1:length(test.x)

    pred = cc_predict(test.D{trl},classifier,1,1);
    x = test.x{trl};

    X = [];
    for st = 1:6
        X(:,st) = smooth(pred(:,st),smoothing);
    end
    
    t1 = X ./ sum(X,2);
    t1(X < (1/6)) = 0;
    
    t2 = t1(find(lags==L)+1:end,:);
    t3 = t1(find(lags==L*2)+1:end,:);
    
    t1 = t1(1:size(t3,1),:);
    t2 = t2(1:size(t3,1),:);    
    
    path1 = [t1(:,1) t2(:,2) t3(:,3)];
    path1(any(path1==0,2),:) = NaN;
    path1 = mean(path1,2);
    
    path1zero = mean([t1(:,1) t2(:,2) t3(:,3)],2);
    
    path2 = [t1(:,4) t2(:,5) t3(:,6)];
    path2(any(path2==0,2),:) = NaN;
    path2 = mean(path2,2);
    
    path2zero = mean([t1(:,4) t2(:,5) t3(:,6)],2);
    
    evidence = [];
    evidence(:,1) = path1 - path2zero;
    evidence(:,2) = path2 - path1zero;
    
    all_evidence(trl,1) = max(evidence(:,1));
    all_evidence(trl,2) = max(evidence(:,2));
        
end

all_evidence(isnan(all_evidence(:,1)),1) = 0;
all_evidence(isnan(all_evidence(:,2)),2) = 0;

[ev,best_trials] = sort(min(all_evidence')','descend');




trl = best_trials(1);

pred = cc_predict(test.D{trl},classifier,1,1);
x = test.x{trl};

X = [];
for st = 1:6
    X(:,st) = smooth(pred(:,st),smoothing);
end

t1 = X ./ sum(X,2);
t1(X < (1/6)) = 0;

t2 = t1(find(lags==L)+1:end,:);
t3 = t1(find(lags==L*2)+1:end,:);

t1 = t1(1:size(t3,1),:);
t2 = t2(1:size(t3,1),:);    

path1 = [t1(:,1) t2(:,2) t3(:,3)];
path1(any(path1==0,2),:) = NaN;
path1 = mean(path1,2);

path1zero = mean([t1(:,1) t2(:,2) t3(:,3)],2);

path2 = [t1(:,4) t2(:,5) t3(:,6)];
path2(any(path2==0,2),:) = NaN;
path2 = mean(path2,2);

path2zero = mean([t1(:,4) t2(:,5) t3(:,6)],2);

evidence = [];
evidence(:,1) = path1 - path2zero;
evidence(:,2) = path2 - path1zero;


figure
sgtitle(['Trial ' num2str(trl)])

subplot(2,1,1)
scatter(1:size(evidence,1),evidence(:,1),'filled'); hold on
scatter(1:size(evidence,1),evidence(:,2),'filled')
title(['Lag = ' num2str(L) ' ms'])

subplot(2,1,2)
imagesc(X')
set(gca,'ticklength',[0 0])
colormap('hot')
ax = gca;
caxis([(1/6) ax.CLim(end)])

