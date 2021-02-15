clear all
clc

rng('shuffle');

%% Parameters

expType = 'MEG'; % 'behav' or 'MEG'

negatorCombos = [
    1 1
    1 2
    1 3
    2 1
    2 2
    2 3
    3 1
    3 2
    3 3
    ];
probabilities = [
    0.1
    0.3
    0.5
    0.7
    0.9
    ];

nptrials = 12;
npblocks = 1;
npforced = 6;

ntrials = 18;
if ~isempty(strmatch(expType,'behav'))
    nblocks = 6;
elseif ~isempty(strmatch(expType,'MEG'))
    nblocks = 10;
end
nforced = 4;

valueRange = [-5:-1, 1:5];
minPathVal = 1;
maxPathVal = 15;

maxNegRow = 3; % can't have more than this many trials in a row with the same negators

safeEV = 1;

nPaths = 2;
nStates = 3;

outcomeConsistency = 0.85; % proportion of the time that one path is better than the other

%% Generate all possible values

nValues = length(valueRange)^(nPaths*nStates); % total number of all possible value combinations

cc = 0;
dcc = round(linspace(1,nValues,100));

values = zeros(nValues,nPaths,nStates); % all possible state values
negValues = zeros(nValues,nPaths,nStates,size(negatorCombos,1)); % all possible negated values (cumulative), per negator combo
EVs = zeros(nValues,size(negatorCombos,1),length(probabilities)); % all possible expected values, per negator combo & probability

t0 = datetime('now');
for p1s1 = 1:length(valueRange)
    for p1s2 = 1:length(valueRange)
        for p1s3 = 1:length(valueRange)
            for p2s1 = 1:length(valueRange)
                for p2s2 = 1:length(valueRange)
                    for p2s3 = 1:length(valueRange)

                        cc = cc + 1;
                        if any(cc == dcc)
                            disp([num2str(round(cc/nValues*100,2)) '% ...']);
                            disp(datetime('now') - t0);
                        end

                        V = [
                            [
                            valueRange(p1s1),...
                            valueRange(p1s2),...
                            valueRange(p1s3)
                            ];...
                            [
                            valueRange(p2s1),...
                            valueRange(p2s2),...
                            valueRange(p2s3)
                            ]
                            ];

                        values(cc,:,:) = V;

                        for n = 1:size(negatorCombos,1)
                            N = negatorCombos(n,:);
                            cV = cumsum(V,2);
                            nV = cV;
                            for path = 1:size(V,1)
                                x = cV(path,1:N(path));
                                if mod(x(end),2) == 1 % if odd
                                    nV(path,N(path)) = nV(path,N(path)) * (-1);
                                    nV(path,N(path):end) = cumsum([nV(path,N(path)) V(path,N(path)+1:end)]);
                                end
                            end
                            negValues(cc,:,:,n) = nV;

                            for p = 1:length(probabilities)
                                EVs(cc,n,p) = sum(nV(:,end)' .* [probabilities(p), 1-probabilities(p)]);
                            end
                        end
                    end
                end
            end
        end
    end
end


% Remove any value combinations that have all even numbers along a path
idx = zeros(nValues,1);
t0 = datetime('now');
for i = 1:nValues
    if any(i == dcc)
        disp([num2str(round(i/nValues*100,2)) '% ...']);
        disp(datetime('now') - t0);
    end
    V = squeeze(values(i,:,:));
    if any(sum(mod(V,2),2) == 0)
        idx(i,1) = 1;
    end
end

values = values(idx == 0,:,:);
negValues = negValues(idx == 0,:,:,:);
EVs = EVs(idx == 0,:,:);

nValues = sum(idx == 0);
disp(['Removed ' num2str(sum(idx)) ' value combinations with all even paths'])


% Remove any value combinations that have a final value outside the range (and excluding 0)
idx = zeros(nValues,1);
for i = 1:nValues
    pathOutcome = squeeze(negValues(i,:,end,:)); % outcome of each path (rows) per negator combo (cols)
    if any(abs(pathOutcome(:)) < minPathVal) || any(abs(pathOutcome(:)) > maxPathVal)
        idx(i,1) = 1;
    end
end

values = values(idx == 0,:,:);
negValues = negValues(idx == 0,:,:,:);
EVs = EVs(idx == 0,:,:);

nValues = sum(idx == 0);
disp(['Removed ' num2str(sum(idx)) ' value combinations with invalid sums'])


% Remove any value combinations that don't give all possible choices, across negators & probabilities
idx = zeros(nValues,1);
minChoice = zeros(nValues,2);
allPathOutcomes = zeros(nValues,2);
uniqueNegators = zeros(nValues,2,2); % value set, path outcome type (aversive 1 or 2), unique negators per path
dcc = linspace(1,nValues,100);
t0 = datetime('now');
for i = 1:nValues

    if any(dcc == i)
        disp([num2str(round(i/nValues,2)*100) '%...'])
        disp(datetime('now') - t0);
    end
    pathOutcome = squeeze(negValues(i,:,end,:)); % outcome of each path (rows) per negator combo (cols)

    % change which path is aversive
    validChoice = zeros(nPaths,2);
    for path = 1:nPaths

        % First see whether there are valid sums for each path pair
        if path == 1
            outcomePerN = pathOutcome(1,:) < 0 & pathOutcome(1,:) < pathOutcome(2,:) & pathOutcome(2,:) > 0; % path 1 is aversive
        elseif path == 2
            outcomePerN = pathOutcome(2,:) < 0 & pathOutcome(2,:) < pathOutcome(1,:) & pathOutcome(1,:) > 0; % path 2 is aversive
        end
        allPathOutcomes(i,path) = sum(outcomePerN);

        if sum(outcomePerN) == 0
            idx(i,1) = 1;
        else
            % check if there is at LEAST 2 unique negators for this path
            uniqueN = negatorCombos(outcomePerN,:);
            uniqueNegators(i,path,:) = [length(unique(uniqueN(:,1))), length(unique(uniqueN(:,2)))];
            % then make sure that a mix of risky and safe decisions can be made at varying probabilities
            for choice = 0:1 % 0 = risky, 1 = safe

                if choice == 0
                    choiceCombos = squeeze(EVs(i,outcomePerN,:)) > safeEV;
                elseif choice == 1
                    choiceCombos = squeeze(EVs(i,outcomePerN,:)) < safeEV;
                end

                validChoice(path,choice+1) = sum(choiceCombos(:)); % how many probability conditions are valid for this choice type

            end
        end
    end
    minChoice(i,path) = min(validChoice(path,:)); % minimum no. of valid probability conditions per choice

    if any(validChoice(:) == 0) % if either choice is impossible, delete
        idx(i,1) = 1;
    end
end

values = values(idx == 0,:,:);
negValues = negValues(idx == 0,:,:,:);
EVs = EVs(idx == 0,:,:);
minChoice = minChoice(idx == 0,:);
allPathOutcomes = allPathOutcomes(idx == 0,:);
uniqueNegators = uniqueNegators(idx == 0,:,:);

nValues = sum(idx == 0);
disp(['Removed ' num2str(sum(idx)) ' value combinations with not enough negator options'])

%% Generate trial structure

iterations = 10;
startingAversive = 2; % 1 or 2

nTrials = (npblocks*nptrials) + (nblocks*ntrials);
tableNames = {'Practice','Block','Trial','Catch','Forced','N','P','V','nV','EV','Optimal'};

% start with most flexible value set
flexiv = find(allPathOutcomes(:,startingAversive) == max(allPathOutcomes(:,startingAversive))); % maximum negator combos that give path 1 (or whatever startingAversive is) as aversive path
flexiv = flexiv(minChoice(flexiv,startingAversive) == max(minChoice(flexiv,startingAversive))); % within the above, find the maximum choice flexibility across probability conditions

ids = randsample(flexiv,iterations,false); % randomly sample no. of iterations from these

allT = cell(iterations,1);

t0 = datetime('now');
dcc = 0;
for it = 1:iterations

    T = array2table(cell(nTrials,length(tableNames)),'VariableNames',tableNames);

    pidx = repmat(probabilities',size(negatorCombos,1),1);
    nidx = repmat([1:size(negatorCombos,1)]',1,length(probabilities));

    % choose starting V (starting with highest number of choice options)
    id = ids(it);
    usedIDs = id;

    % practice/block/trial info
    T.Practice = [ones(npblocks*nptrials,1); zeros(nblocks*ntrials,1)];
    x1 = repmat(1:npblocks,nptrials,1);
    x2 = repmat(1:nblocks,ntrials,1);
    T.Block = [x1(:); x2(:)];
    T.Trial = [repmat(1:nptrials,1,npblocks)'; repmat(1:ntrials,1,nblocks)'];

    for cc = 1:size(T,1)

        dcc = dcc + 1;
        disp(['Iteration #' num2str(it) ', ' num2str(round((dcc-1)/(iterations*nTrials))*100) '% ...']);
        disp(datetime('now') - t0);

        % which path (1 or 2) is the aversive path for this block
        if T.Practice(cc) == 0 && T.Block(cc) > max(T.Block)/2
            aversivePath = setdiff(1:2,startingAversive);
        else
            aversivePath = startingAversive;
        end
        rewardingPath = setdiff([1 2],aversivePath);

        % whether this is a forced trial or not
        if T.Trial(cc) == 1
            if T.Practice(cc) == 1
                forcedPool = randsample([ones(npforced/2,1);ones(npforced/2,1)*2],npforced,false);
            else
                forcedPool = randsample([ones(nforced/2,1);ones(nforced/2,1)*2],nforced,false);
            end
        end

        if (T.Practice(cc) == 1 && T.Trial(cc) <= npforced) || (T.Practice(cc) == 0 && T.Trial(cc) <= nforced)
            forcedPool
            thisForced = forcedPool(randsample(1:length(forcedPool),1));
            thisForced
            T.Forced(cc) = {thisForced};
            forcedPool(randsample(find(forcedPool == thisForced),1),:) = [];
        else
            T.Forced(cc) = {0};
        end

        % which values for this block
        if cc == 1
            V = squeeze(values(id,:,:));
        elseif T.Trial(cc) == 1 % change the values at the start of each block

            disp(['........ calculating new value set for block ' num2str(T.Block(cc)) '...']);

            % find value set with one different value in each path
            vsimilarity = zeros(size(values,1),1);
            for v = 1:size(values,1)
                compv = sum(V == squeeze(values(v,:,:)),2);
                if compv(1,1) == 2 && compv(2,1) == 2
                    vsimilarity(v,1) = 1;
                end
            end

            for i = 1:length(usedIDs)
                vsimilarity(usedIDs(i),1) = 0; % don't include value sets that have already been used
            end

            % narrow down to sets with the most different negator combos
            vsimnegators = zeros(nValues,size(negatorCombos,1));
            for a = 1:nValues
                pathOutcome = squeeze(negValues(a,:,end,:));
                vsimnegators(a,:) = pathOutcome(aversivePath,:) < pathOutcome(rewardingPath,:) & ...
                                pathOutcome(aversivePath,:) < 0 & ...
                                pathOutcome(rewardingPath,:) > 0;
            end

            tmpn = cell2mat(T.N);
            negatorsSoFar = zeros(size(negatorCombos,1),1);
            for n = 1:length(negatorsSoFar)
                negatorsSoFar(n,1) = sum(tmpn(:,1) == negatorCombos(n,1) & tmpn(:,2) == negatorCombos(n,2)) / size(tmpn,1);
            end
            rareNegators = 1 - negatorsSoFar;

            validv = [sum(vsimnegators,2) sum(vsimnegators .* rareNegators',2)];
            [validv,validvidx] = sortrows(validv(vsimilarity == 1,:),'descend');

            validvidx = validvidx(validv(:,1) == validv(1,1) & validv(:,2) == validv(1,2),:);
            vsimilarityidx = find(vsimilarity);
            validvidx = vsimilarityidx(validvidx);

            % randomly select and update V
            id = randsample(validvidx,1);
            usedIDs = [usedIDs, id];

            V = squeeze(values(id,:,:));
        end

        if T.Trial(cc) > 1
            lastN = cell2mat(T.N(T.Practice == T.Practice(cc) & T.Block == T.Block(cc) & T.Trial < T.Trial(cc)));
            lastP = cell2mat(T.P(T.Practice == T.Practice(cc) & T.Block == T.Block(cc) & T.Trial < T.Trial(cc)));
        end

        % will this trial be normal (i.e. one path aversive, one path rewarding) or just random (i.e. catch trial)?
        if T.Trial(cc) == 1
            catchTrial = false;
        elseif cc > 1
            if cell2mat(T.Catch(cc-1)) ~= 1 % can't have 2 catch trials in a row
                catchTrial = randsample([false true],1,true,[outcomeConsistency 1-outcomeConsistency]);
            else
                catchTrial = false;
            end
        end

        if ~catchTrial
            pathOutcome = squeeze(negValues(id,:,end,:));
            validAversive = pathOutcome(aversivePath,:) < pathOutcome(rewardingPath,:) & ...
                pathOutcome(aversivePath,:) < 0 & ...
                pathOutcome(rewardingPath,:) > 0;
        else
            validAversive = 1:size(pidx,1);
        end

        this_pidx = pidx(validAversive,:);
        this_nidx = nidx(validAversive,:);

        if cell2mat(T.Forced(cc)) == 1
            fidx = probabilities >= 0.5;
        elseif cell2mat(T.Forced(cc)) == 2
            fidx = probabilities <= 0.5;
        else
            fidx = ones(length(probabilities),1);
        end

        % randomly choose correct choice for this trial
        choice = randsample([0 1],1); % 0 = risk, 1 = safe

        if cell2mat(T.Forced(cc)) == 0
            if choice == 0
                choiceCombos = squeeze(EVs(id,validAversive,:)) > safeEV;
            elseif choice == 1
                choiceCombos = squeeze(EVs(id,validAversive,:)) < safeEV;
            end
        else
            rc = squeeze(EVs(id,validAversive,fidx)) > safeEV;
            lc = squeeze(EVs(id,validAversive,fidx)) < safeEV;
            if sum(rc(:)) == 0 && sum(lc(:)) > 0
                choice = 1;
                choiceCombos = lc;
            elseif sum(rc(:)) > 0 && sum(lc(:)) == 0
                choice = 0;
                choiceCombos = rc;
            elseif sum(rc(:)) > 0 && sum(lc(:)) > 0
                if randsample([0 1],1) == 0
                    choice = 0;
                    choiceCombos = rc;
                else
                    choice = 1;
                    choiceCombos = lc;
                end
            else
                error('Invalid choice combinations for forced choice trial!');
            end
        end

        % randomly choose from negator-probability array (but try to avoid previous selections)
        if T.Trial(cc) > 1

            % weight probability choice by experiment history
            choosep = ones(length(probabilities),1)*length(lastP);
            for p = 1:length(probabilities)
                choosep(p,1) = choosep(p,1) - sum(lastP == probabilities(p));
            end
            if cell2mat(T.Forced(cc)) > 0
                choosep(fidx == 0,:) = 0;
            end
            choosep = choosep / max(choosep);
            choosep = choosep / sum(choosep);

            % weight negator choice probability by experiment history
            choosen = ones(size(negatorCombos,1),1)*size(lastN,1);
            choosen(setdiff(1:size(negatorCombos,1),unique(this_nidx)),1) = 0; % remove invalid negator combos
            for n = 1:size(negatorCombos,1)
                choosen(n,1) = choosen(n,1) - sum(lastN(:,1) == negatorCombos(n,1) | lastN(:,2) == negatorCombos(n,2))/2; % down-weight if one of the negators has been used for that path already
                choosen(n,1) = choosen(n,1) - sum(lastN(:,1) == negatorCombos(n,1) & lastN(:,2) == negatorCombos(n,2))/2; % down-weight if exact combo has been used
            end
            choosen(setdiff(1:size(negatorCombos,1),unique(this_nidx)),:) = []; % remove invalid negator combos again
            choosen = choosen / max(choosen);
            choosen = choosen / sum(choosen);

            % combine weights
            cw = repmat(choosep',length(unique(this_nidx)),1);
            cw  = cw + repmat(choosen,1,length(probabilities));
            cw = cw ./ sum(cw(:));
            cw(choiceCombos == 0) = 0;

            % select
            npid = randsample(prod(size(cw)),1,'true',cw(:));

        else
            npid = randsample(find(choiceCombos),1);
        end

        N = negatorCombos(this_nidx(npid),:);
        P = this_pidx(npid);
        nV = squeeze(negValues(id,:,:,this_nidx(npid)));
        EV = squeeze(EVs(id,this_nidx(npid),find(probabilities == P)));

        % fill in table
        T.N(cc) = {N};
        T.P(cc) = {P};
        T.V(cc) = {V};
        T.nV(cc) = {nV};
        T.EV(cc) = {EV};
        T.Catch(cc) = {catchTrial};
        if EV > safeEV
            T.Optimal(cc) = {'risk'};
        elseif EV < safeEV
            T.Optimal(cc) = {'safe'};
        elseif EV == safeEV
            T.Optimal(cc) = {'either'};
        end

    end

    % Add balanced transition outcomes to forced choice trials
    for block = 1:npblocks % practice
        idx1 = find(T.Practice == 1 & T.Block == block);
        idx2 = find(cell2mat(T.Forced(idx1)));
        F = randsample([ones(length(idx2)/2,1); ones(length(idx2)/2,1)*2],length(idx2),false);
        for f = 1:length(F)
            T.Forced(idx1(idx2(f))) = {F(f)};
        end
    end
    for block = 1:nblocks % task
        idx1 = find(T.Practice == 0 & T.Block == block);
        idx2 = find(cell2mat(T.Forced(idx1)));
        F = randsample([ones(length(idx2)/2,1); ones(length(idx2)/2,1)*2],length(idx2),false);
        for f = 1:length(F)
            T.Forced(idx1(idx2(f))) = {F(f)};
        end
    end

    allT{it} = T;

end

% show negator frequencies across iterations
evenness = zeros(iterations,2);
otherCounts = zeros(iterations,2); % catch trial rate, average EV

cmap = [141 51 255; 51 221 255; 255 162 51]/255;
figure
cc = 0;
for it = 1:iterations

    cc = cc + 1;
    subplot(2,iterations,cc);
    pdata = cell2mat(allT{it}.N);
    imagesc(pdata);
    title(['Iteration ' num2str(it)]);
    ylabel('trials')
    xlabel('paths')
    colormap(cmap);

    evenness(it,:) = sum(diff(pdata) ~= 0);

    cc = cc + 1;
    subplot(2,iterations,cc);
    pdata = cell2mat(allT{it}.N);
    pdata = [sum(pdata(:,1) == 1), sum(pdata(:,1) == 2), sum(pdata(:,1) == 3);...
        sum(pdata(:,2) == 1), sum(pdata(:,2) == 2), sum(pdata(:,2) == 3)] / size(pdata,1);
    B = bar(pdata);
    for b = 1:length(B)
        B(b).FaceColor = cmap(b,:);
    end
    hold on
    ax = gca;
    plot(ax.XLim,[.5 .5],'k--');
    ylim([0 1])
    xlabel('path')
    ylabel('proportion')

    otherCounts(it,1) = mean(cell2mat(allT{it}.Catch));
    otherCounts(it,2) = mean(cell2mat(allT{it}.EV));

end

[~,bestIteration] = sort(mean(evenness,2),'descend');

% check proportion of safe/risky choices
choicePcnt = zeros(iterations,1);
for it = 1:iterations
    tmp = allT{1}.Optimal;
    choicePcnt(it,1) = length(find(not(cellfun('isempty',strfind(tmp,'risk')))))/length(tmp);
end

%% Plot

close all
it = 1;
T = allT{it};

% structure
cmap = [141 51 255; 51 221 255; 255 162 51]/255;
figure
imagesc(cell2mat(allT{it}.N));
title(['Iteration ' num2str(it)]);
ylabel('trials')
xlabel('paths')
colormap(cmap);

figure
pdata = cell2mat(allT{it}.N);
pdata = [sum(pdata(:,1) == 1), sum(pdata(:,1) == 2), sum(pdata(:,1) == 3);...
    sum(pdata(:,2) == 1), sum(pdata(:,2) == 2), sum(pdata(:,2) == 3)] / size(pdata,1);
B = bar(pdata);
for b = 1:length(B)
    B(b).FaceColor = cmap(b,:);
end
hold on
ax = gca;
plot(ax.XLim,[.5 .5],'k--');
ylim([0 1])
xlabel('path')
ylabel('proportion')

% negators
pdata = cell2mat(T.N);
figure
for i = 1:nPaths
    title(['Negators for path ' num2str(i)])
    subplot(1,nPaths,i);
    plot(pdata(:,i)); hold on
    scatter(1:size(pdata,1), pdata(:,i), 'filled');
end
xlim([1 size(pdata,1)])

% probabilities
pdata = cell2mat(T.P);
F = figure
title('Probabilities')
plot(pdata,'Color',[1 1 1]); hold on
cmap = [144 12 63; 199 0 57; 255 87 51; 255 195 0; 255 255 50]/255;
for p = 1:length(probabilities)
    scatter(find(pdata == probabilities(p)), ones(sum(pdata == probabilities(p)),1)*probabilities(p),...
        'filled', 'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(p,:));
end
darkBackground(F,[0 0 0],[1 1 1]);
xlim([1 size(pdata,1)])
set(gca,'TickLength',[0 0])

% ev
pdata = cell2mat(T.EV);
figure
subplot(2,1,1)
title('EV')
plot(pdata); hold on
scatter(1:size(pdata,1), pdata, 'filled');
plot([1 size(pdata,1)],[safeEV safeEV],'k');
xlim([1 size(pdata,1)])

subplot(2,1,2) % choice proportion
bar([sum(pdata > safeEV), sum(pdata == safeEV), sum(pdata < safeEV)] / length(pdata)); hold on
plot([0 4],[.5 .5],'k')
set(gca,'XTickLabels',{'Risky','Either','Safe'})
xlabel('Optimal choices')
ylabel('Percentage')
xlim([0 4])
ylim([0 1])

% state values
pdata = zeros(nPaths,nStates,length(usedIDs));
for i = 1:length(usedIDs)
    pdata(:,:,i) = squeeze(values(usedIDs(i),:,:));
end

figure
for path = 1:nPaths
    for st = 1:nStates
        plot(1:length(usedIDs),squeeze(pdata(path,st,:))); hold on
        scatter(1:length(usedIDs),squeeze(pdata(path,st,:)),'filled');
    end
end
title('State values')
ylabel('value')
xlabel('block')
ylim([valueRange(1)-1 valueRange(end)+1])

% end values
pdata = zeros(size(T,1),2);
for trl = 1:size(T,1)
    pdata(trl,1) = T.nV{trl}(1,end);
    pdata(trl,2) = T.nV{trl}(2,end);
end

figure
plot(pdata);
hold on
ax = gca;
plot(ax.XLim,[safeEV safeEV])
title('Path values')

%% Write JSON file

apostrophe = '''';

json = strcat('var exp_list = ',apostrophe,'[{');
for col = 1:size(T,2)
    json = strcat(json,'"',tableNames{col},'":{');
    for row = 0:size(T,1)-1
        x = [];
        if col < size(T,2)
            if iscell(table2array(T(row+1,col)))
                x = cell2mat(table2array(T(row+1,col)));
            else
                x = table2array(T(row+1,col));
            end
            if length(x) == 1
                x = num2str(x);
            elseif size(x,1) == 1 && size(x,2) > 1
                y = '[';
                for i = 1:size(x,2)
                    if i < size(x,2)
                        y = strcat(y,num2str(x(i)),',');
                    else
                        y = strcat(y,num2str(x(i)),']');
                    end
                end
                x = y;
            elseif size(x,1) > 1 && size(x,2) > 1
                y = '[';
                for i = 1:size(x,1)
                    y = strcat(y,'[');
                    for j = 1:size(x,2)
                        if j < size(x,2)
                            y = strcat(y,num2str(x(i,j)),',');
                        else
                            if i == 1
                                y = strcat(y,num2str(x(i,j)),'],');
                            else
                                y = strcat(y,num2str(x(i,j)),']]');
                            end
                        end
                    end
                end
                x = y;
            end
        else
            if strmatch('safe',table2array(T(row+1,col)))
                x = num2str(1);
            elseif strmatch('risk',table2array(T(row+1,col)))
                x = num2str(0);
            end
        end
        json = strcat(json,'"',num2str(row),'":',x);
        if row < size(T,1)-1
            json = strcat(json,',');
        elseif row == size(T,1)-1 && col < size(T,2)
            json = strcat(json,'},');
        else
            json = strcat(json,'}}]',apostrophe);
        end
    end
end

% Y = {};
% cc = 0;
% for i = 1:999:length(json)
%     cc = cc + 1;
%     try
%         Y{cc} = json(i:i+998);
%     catch
%         Y{cc} = json(i:end);
%     end
% end
% Y = Y';

if startingAversive == 1
    label = [expType '-A'];
elseif startingAversive == 2
    label = [expType '-B'];
end

writetable(T,['structure_' label '.csv']);
tmp = readtable(['structure_' label '.csv']); % swap the columns back around the right way
tmp(:,9:20) = array2table(table2array([tmp(:,9:2:14), tmp(:,10:2:14), tmp(:,15:2:20), tmp(:,16:2:20)]));
writetable(tmp,['structure_' label '.csv']);

fid = fopen(['structure_' label '.js'],'wt');
fprintf('%s',json);
fclose(fid);
    
