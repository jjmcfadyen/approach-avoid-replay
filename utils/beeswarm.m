function [x,y] = beeswarm(y,binwidth,width)
% y = vector of data
% binwidth = resolution at which to bin the data into swarm
% width = width of spread along x axis (0-1)

if width < 0 || width > 1
    error('Width must be between 0 and 1')
end

N = length(y);

origy = y;
fullx = zeros(N,1);

[y,sortidx] = sort(origy);
x = zeros(N,1);

% see which data points have nearby points
neighbours = zeros(N,N);
for i = 1:N
    idx = setdiff(find(y > y(i)-(binwidth/2) & y < y(i)+(binwidth/2)),i);
    neighbours(i,idx) = 1;
end
num_neighbours = sum(neighbours,2);

% only look at data that has neighbours
pointidx = num_neighbours > 0;
y = y(pointidx);
x = x(pointidx);

% see how many data points in each bin
steps = NaN;
thisy = y;
factorincrease = 0;
while length(steps)==1
    steps = linspace(min(thisy),max(thisy),ceil((max(thisy)-min(thisy))/binwidth));
    thisy = thisy*10;
    factorincrease = factorincrease + 1;
end
steps = steps * 10^(factorincrease-1);

n = nan(length(steps)-1,1);
for i = 1:length(steps)-1
    idx = find(y >= steps(i) & y < steps(i+1));
    n(i,1) = length(idx);
end

% assign to x axis
maxn = max(n);
if maxn/2 - round(maxn/2) == 0 % even
    even_dist = linspace(-width,width,maxn);
    odd_dist = linspace(-width,width,maxn-1);
else
    even_dist = linspace(-width,width,maxn-1);
    odd_dist = linspace(-width,width,maxn);
end
for i = 1:length(n)
    idx = find(y >= steps(i) & y < steps(i+1));
    if n(i) > 1
        if n(i)/2 - round(n(i)/2) == 0 % even
            startidx = length(even_dist)/2-(n(i)/2)+1;
            x(idx,:) = even_dist(startidx:(startidx+n(i)-1));
        else % odd
            startidx = ceil(length(odd_dist)/2)-floor(n(i)/2);
            x(idx,:) = odd_dist(startidx:(startidx+n(i)-1));
        end
        x(idx,:) = randsample(x(idx,:),length(idx),false);
    end
end

% put back into x and y
fullx(pointidx) = x;
x = fullx(sortidx,:);
y = origy;

end