meanamp = [];
cidx = [];
for trl = 1:nTrls
    
    cfg = [];
    cfg.trials = trl;
    thistrial = ft_selectdata(cfg,data); % 100 hz data for replay
    
    [onsets, seqevidence] = get_replayOnsets(thistrial,classifier,lagrange);
    
    seqevidence = seqevidence{1}(sum(isnan(seqevidence{1}),2)==0,:);
    path1 = mean(seqevidence(:,1:2),2);
    path2 = mean(seqevidence(:,3:4),2);
    x = thistrial.time{1}(1:size(seqevidence,1));
    
    xidx = x>=0 & x<=5;
    path1 = path1(xidx);
    path2 = path2(xidx);
    x = x(xidx);
    
    if behav.Forced(trl)==0 && ((behav.nV_1(trl)>0 && behav.nV_2(trl)<0) || (behav.nV_1(trl)<0 && behav.nV_2(trl)>0)) && ~any(isnan(x))
        if behav.nV_1(trl) > behav.nV_2(trl)
            meanamp(trl,1,1:length(x)) = path1;
            meanamp(trl,2,1:length(x)) = path2;
        else
            meanamp(trl,1,1:length(x)) = path2;
            meanamp(trl,2,1:length(x)) = path1;
        end
        cidx = [cidx; behav.Choice(trl)];
    end
end

figure
for c = 1:2
    subplot(1,2,c)
    for p = 1:2
        m = squeeze(mean(meanamp(cidx==c,p,:)))';
        smoothm = smooth(m,20);
%         sem = squeeze(std(meanamp(cidx==c,p,:))/sqrt(size(meanamp,1)))';
%         upper = m+sem;
%         lower = m-sem;
%         patch([x fliplr(x)],[upper fliplr(lower)],cmap(p,:),'facealpha',.2,'edgealpha',0); hold on
        plot(x,m,'linewidth',1.3,'color',cmap(p,:),'linestyle',':'); hold on
        plot(x,smoothm,'linewidth',1.6,'color',cmap(p,:)); hold on
    end
    xlim([1 5])
    set(gca,'ticklength',[0 0])
end



