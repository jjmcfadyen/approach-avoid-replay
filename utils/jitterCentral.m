function xpos = jitterCentral(x,b,w,i)
% x = data (vector)
% b = no. of bins (for histogram)
% w = width
% i = centre point

n = length(x);

[N,edges] = histcounts(x,b);

xpos = nan(n,1);
for e = 1:length(edges)
    if e < length(edges)
        L = edges([e e+1]);
    else
        L = [edges(e) max(x)];
    end
    this_x = find(x > L(1) & x < L(2));
    if length(this_x) == 1
        xpos(this_x) = i;
    else
        xpos(this_x) = linspace(i-(w/2),i+(w/2),length(this_x));
    end
end

end