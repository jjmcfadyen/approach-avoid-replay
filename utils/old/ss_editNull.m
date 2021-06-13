function I = ss_editNull(U,opts)

% U = list of null permutations (including initial non-null permutation)
% I = index of permutations to keep

if ~isfield(opts,'ignoreOrder')
    opts.ignoreOrder = false;
end
if ~isfield(opts,'removePos')
    opts.removePos = false;
end
if ~isfield(opts,'removeLoop')
    opts.removeLoop = false;
end

states = unique(U);
I = [1:size(U,1)]';

%% Ignore order
% If looking at general sequenceness, then all transitions are averaged
% together, meaning that the order of the transitions does not matter

if opts.ignoreOrder
    uTrans = [nchoosek(states,2); fliplr(nchoosek(states,2))];
    uIdx = zeros(size(U,1),size(uTrans,1));
    for perm = 1:length(I)
        for t = [1 2 4 5]
            uIdx(perm,ismember(uTrans,U(I(perm),[t t+1]),'rows')) = 1;
        end
    end

    [~,I,~] = unique(uIdx,'rows');
    I = sort(I);
end

%% Remove transitions to same position across paths
% e.g. state 1 to 4
% e.g. state 2 to 5
% e.g. state 3 to 6

if opts.removePos
    ridx = [];
    for perm = 1:length(I)
        for t = [1 2 4 5]
            thisTrans = U(I(perm),[t t+1]);
            if abs(diff(thisTrans)) == 3
                ridx = [ridx; perm];
            end
        end
    end
    I(unique(ridx),:) = [];
end

%% Remove transitions that loop back to the other path
% e.g. state 6 to 1
% e.g. state 1 to 6

if opts.removeLoop
    ridx = [];
    for perm = 1:length(I)
        for t = [1 2 4 5]
            thisTrans = sort(U(I(perm),[t t+1]));
            if ismember(thisTrans,[1 6],'rows')
                ridx = [ridx; perm];
            end
        end
    end
    I(unique(ridx),:) = [];
end

end