function [U] = generate_nullperms(opts)
% Generate the "null" transitions for permutation testing

if nargin == 1
    if ~isfield(opts,'pathLink') % 3 to 4 transition
        opts.pathLink = false;
    end
    if ~isfield(opts,'withinSeq') % 1 to 3 or 4 to 6
        opts.withinSeq = false;
    end
    if ~isfield(opts,'removeMirror') % if agnostic about the WITHIN-path transitions, remove duplicates
        opts.removeMirror = false;
    end
end

%% Parameters

nStates = 6;

% Outline all possible transitions
transitions = [1 2; 2 3; 4 5; 5 6];

% Null transitions (all cross-sequence, except 3-4 pair)
nt = [
    1 4
    1 5
    1 6
    2 4
    2 5
    2 6
    3 5
    3 6
    ];

if opts.pathLink
    nt = [nt; 3 4];
end
if opts.withinSeq
    nt = [nt; 1 3; 4 6];
end

% Get all possible orders of states
[~, U] = uperms(1:nStates);
U = U(2:end,:); % remove first row (non-null transitions: 1 to nStates)

% Get rid of any permutation that has a non-null transition
ridx = zeros(size(U,1),1); % permutations to remove
for perm = 1:size(U,1)
    thisPerm = U(perm,:);
    for t = 1:size(transitions,1)
        thisTrans = thisPerm(transitions(t,1:2));
        validF = find(nt(:,1) == thisTrans(1) & nt(:,2) == thisTrans(2)); % forward
        validB = find(nt(:,2) == thisTrans(1) & nt(:,1) == thisTrans(2)); % backward
        if isempty(validF) && isempty(validB)
            ridx(perm,1) = 1;
        end
    end
end
U = U(ridx == 0,:);

% Remove mirrored paths
if opts.removeMirror
    
    TM = nan(size(U,1),size(transitions,1),size(transitions,2));
    for perm = 1:size(U,1)
        thisPerm = U(perm,:);
        for t = 1:size(transitions,1)
            TM(perm,t,:) = thisPerm([transitions(t,1),transitions(t,2)]);
        end
    end
    
    duplicates = [];
    for perm1 = 1:size(U,1)
        for perm2 = 1:size(U,1)
            if perm1 ~= perm2
                t1 = squeeze(TM(perm1,:,:));
                t2 = squeeze(TM(perm2,:,:));
                for r = 1:size(t1,1)
                    t1(r,:) = sort(t1(r,:));
                    t2(r,:) = sort(t2(r,:));
                end
                t1 = squash(sortrows(t1));
                t2 = squash(sortrows(t2));
                if all(t1==t2)
                    duplicates = [duplicates; sort([perm1 perm2])];
                end
            end
        end
    end
    duplicates = unique(duplicates,'rows');
    
    U = U(unique(duplicates(:,1)),:);
    
end

% Add actual transitions back in
U = [1:nStates; sortrows(U)];

end