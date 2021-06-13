function [M, thresh] = so_subtractNull(x,g,stype)

% x = iterations (trials, subjects) x permutations x lags

if length(size(x))==2 % no iterations (trials, subjects)
    oldx = x;
    x = [];
    x(1,:,:) = oldx;
    stype = 'trial';
end

switch stype
    case 'trial'
        nTrls = size(x,1);
    case 'group'
        nTrls = 1;
end

M = nan(nTrls,size(x,3));
thresh = nan(nTrls,1);
for trl = 1:nTrls
    
    switch stype
        case 'trial'
            y = squeeze(x(trl,:,:));
        case 'group'
            y = squeeze(mean(x));
    end
    
    m = y(1,:);
    np = y(2:end,:);
    npThresh = quantile(max(np,[],2),.95);
    
    if g == 3
        m = abs(m);
        np = abs(np);
        npThresh = quantile(max(np,[],2),.975);
    end
    
    M(trl,:) = m - npThresh;
    thresh(trl,:) = npThresh;
    
end

end