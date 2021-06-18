function best_traintime = optimise_replay(replay,criteria_seqType,thesetimes,CV)

% Determine criteria for each training time
sigreplay = [];
for g = 1:3
    for t = 1:size(replay,1)
        opts = [];
        opts.makeplot = false;
        opts.subtractNull = true;
        opts.g = g;
        tmp = ss_plot(squeeze(replay(t,:,:,g,:)),opts);
        sigreplay(t,g,:) = tmp.m;
    end
end

zreplay = squeeze(mean(replay(:,:,1,:,:),2));
mcreplay = squeeze(mean(replay(:,:,1,:,:),2));
for g = 1:3
    for i = 1:size(zreplay,1)
        zreplay(i,g,:) = zscore(zreplay(i,g,:));
        mcreplay(i,g,:) = mcreplay(i,g,:) - mean(mcreplay(i,g,:));
    end
end

criteria = array2table(nan(length(thesetimes),9),'variablenames',{'peakSig','peakZ','peakVal','startVal','cvAccuracy','seqType','peakLocF','peakLocB','peakLocD'});
for t = 1:length(thesetimes)
    
    d = squeeze(mean(replay(t,:,1,:,:),2));
    
    % get either the maximum forward OR backward value
    switch criteria_seqType
        case 'either'
            [pk,loc] = max(squeeze(sigreplay(t,1:2,:)),[],2);
            criteria.peakLocF(t) = loc(1);
            criteria.peakLocB(t) = loc(2);
            
            idx = find(pk==max(pk));
            criteria.seqType(t) = idx;
            
            pk = max(pk);
            loc = loc(pk==max(pk));
        case 'diff'
            [pk,loc] = max(squeeze(sigreplay(t,3,:)));
            criteria.peakLocD(t) = loc;
            idx = 3;
            d = abs(d);
    end
    
    criteria.peakSig(t) = pk;
    criteria.peakVal(t) = squeeze(d(idx,loc));
    criteria.peakZ(t) = abs(zreplay(t,idx,loc));
    criteria.startVal(t) = abs(d(idx,1));
    criteria.cvAccuracy(t) = CV(t);
end

% get median forward and backward replay
switch criteria_seqType
    case 'either'
        medForward = round(median(criteria.peakLocF));
        medBackward = round(median(criteria.peakLocB));
        
        criteria.peakLocDiff = nan(size(criteria,1),1);
        criteria.peakLocDiff(criteria.seqType==1) = abs(criteria.peakLocF(criteria.seqType==1)-medForward);
        criteria.peakLocDiff(criteria.seqType==2) = abs(criteria.peakLocB(criteria.seqType==2)-medBackward);
        
    case 'diff'
        criteria.peakLocDiff = abs(criteria.peakLocD - median(criteria.peakLocD));
end

criteria = criteria(:,[3 4]);
bestCombination = [max(criteria.peakVal) min(criteria.startVal)];

% Pick best match
[~,best_traintime] = min(pdist2(table2array(criteria),bestCombination));

end