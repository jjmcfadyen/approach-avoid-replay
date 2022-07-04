function [startvals,fitvals,err] = cluster_fit_startrange(filename)
% [startvals,fitvals,err] = cluster_fit_startrange(filename)
% fits model using a range of starting vlaues

%% Load file and get details

tmp = load(filename);

model = tmp.model;
d = tmp.thisd;
numStart = tmp.numStart;
info = tmp.info;

[FILEPATH,NAME,EXT] = fileparts(filename);

% if contains(tmp.model.name,'random')
%     nIterations = 100;
% else
    nIterations = 1;
% end

%% Run function

% get function for this model
fun_model = @(x) sim_model(d,model,x,info.randtype);

% get lower and upper bounds per parameter
pnames = fieldnames(model.params);
np = length(pnames);
freeidx = zeros(1,length(pnames));
for p = 1:np
    freeidx(p) = isnan(model.params.(pnames{p}).val);
end

pnames = pnames(freeidx==1);
np = length(pnames);

lb = nan(1,np);
ub = nan(1,np);
for p = 1:np
    lb(p) = model.params.(pnames{p}).lb;
    ub(p) = model.params.(pnames{p}).ub;
end

% extract starting values of free parameters
while numStart^np > 100
    warning(['Too many combinations of starting values to compute (' num2str(numStart^np) '). Reducing numStart down by 1.' ])
    numStart = numStart - 1;
end

startvals = nan(np,numStart);
for p = 1:np
    if numStart==1
        startvals(p,:) = model.params.(pnames{p}).start;
    else
        if isinf(ub(p))
            thisrange = exprnd(4,1,numStart-1);
        else
            thisrange = linspace(lb(p),ub(p),numStart-1);
        end
        startvals(p,:) = sort([model.params.(pnames{p}).start thisrange]);
    end
end

if size(startvals,2)==2
    startvals = combvec(startvals(1,:),startvals(2,:));
elseif size(startvals,2)==3
    startvals = combvec(startvals(1,:),startvals(2,:),startvals(3,:));
elseif size(startvals,2)==4
    startvals = combvec(startvals(1,:),startvals(2,:),startvals(3,:),startvals(4,:));
elseif size(startvals,2)>4
    error(['Too many combinations of starting values to compute (' num2str(numStart^np) ')' ])
end

%% Optimise free parameters for each combination of starting values

fitvals = [];
err = [];
for it = 1:nIterations

    optstart = nan(size(startvals,2),np+1);
    for v = 1:size(startvals,2)
        if numStart>1
            disp(['--- Model ' model.name ': optimising parameters for start set ' num2str(v) ' of ' num2str(size(startvals,2)) '...'])
        end
        [x,err] = patternsearch(fun_model,startvals(:,v),[],[],[],[],lb,ub);
        optstart(v,:) = [err x'];
    end
    
    err(it) = optstart(:,1)';
    fitvals{it} = array2table(optstart(:,2:end),'variablenames',pnames');

end

% Rearrange
if nIterations==1
    fitvals = fitvals{it};
else
    tmp = [];
    for it = 1:nIterations
        tmp(it,:,:) = table2array(fitvals{it});
    end
    fitvals = array2table(squeeze(mean(tmp))','variablenames',fitvals{1}.Properties.VariableNames);
end

% save
save(fullfile(FILEPATH,[NAME EXT]),'startvals','fitvals','err','d','model','info');

end