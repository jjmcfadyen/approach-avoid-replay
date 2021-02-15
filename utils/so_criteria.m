function criteria = so_criteria(allTraces,grand,opts)

nLambda = size(allTraces,1);

% get max values
peakLoc = nan(nLambda,3);
peakVal = nan(nLambda,3);
meanVal = nan(nLambda,3);
for g = 1:3

    % get the grand average mean
    gm = squeeze(mean(grand(:,g,:)));

    % determine the peaks of interest
    [gpks, glocs] = so_findpeaks(gm,g);

    % Get maximum peak & location
    [gmax,gloc] = max(gpks);
    gloc = glocs(gloc,1);
    
    if isempty(gloc)
        gloc = 0;
        gmax = max(gm);
    end

    for i = 1:nLambda

        y = squeeze(mean(allTraces(i,:,g,:)));

        [pk, loc] = so_findpeaks(y,g);
        [ymax,yloc] = max(pk);
        yloc = loc(yloc,1);
        
        if isempty(yloc)
            yloc = 0;
            ymax = max(y);
        end

        peakLoc(i,g) = abs(yloc - gloc);
        peakVal(i,g) = ymax - gmax;
        meanVal(i,g) = max(y) - min(y); % range
        startVal(i,g) = squeeze(mean(allTraces(i,:,1,g,1)));
    end
end

criteria.peakLoc  = peakLoc;
criteria.peakVal  = peakVal;
criteria.meanVal  = meanVal;
criteria.startVal = startVal;

end