function [X,Y] = organise_FL_data(data,traintime,nulldata)
% Organise data from functional localiser for building classifiers

%% Get info from data

nTrls = length(data.trial);
nChan = length(data.label);

Y = data.trialinfo; % labels

%% Organise data

% Subset just the one timepoint from each trial/channel
X = nan(nTrls,nChan);
for trl = 1:nTrls
    X(trl,:) = data.trial{trl}(:,findMin(traintime/1000,data.time{trl})); 
end

% Scale the data matrix
X = X ./ prctile(abs(X(:)),95);

%% Add null data

if nulldata>0
    X = [X; zeros(nTrls,nChan)];
    Y = [Y; zeros(nTrls,1)];
end

end