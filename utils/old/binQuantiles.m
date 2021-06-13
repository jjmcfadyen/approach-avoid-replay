function [y idx edges] = binQuantiles(x,nbin)

if size(x,2) > size(x,1)
    x = x'; % one column of many rows
end

I = size(x,1);

q = quantile(x,linspace(0,1,nbin+1));

edges = cell(nbin,1);
for i = 1:nbin
    edges{i,1} = [num2str(q(i)) ' to ' num2str(q(i+1))];
end

Q = q;
q(1) = -Inf;

y = nan(size(x,1),1);
idx = nan(size(x,1),1);
for i = 1:I
    idx(i,1) = find(x(i,1) > q,1,'last');
    y(i,1) = Q(idx(i,1));
end

end