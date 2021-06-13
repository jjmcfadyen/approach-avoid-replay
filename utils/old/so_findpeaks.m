function [P, L] = so_findpeaks(y,g)

w = 3; % minimum width of peaks

% Smooth y
sy = smooth(y);

% First get positive peaks
[P1,L1,W] = findpeaks(sy);
idx = round(W) >= w;
P1 = P1(idx); % remove tiny peaks (smaller than 'w' width parameter)
L1 = L1(idx);

% Then negative peaks
[P2,L2,W] = findpeaks(-sy);
idx = round(W) >= w;
P2 = P2(idx); % remove tiny peaks (smaller than 'w' width parameter)
L2 = L2(idx);

% Output
if g < 3
    P = y(L1); % positive peaks
    L = L1;
else
    [L,sortidx] = sort([L1; L2]); % both peaks
    P = y(L);
end

end