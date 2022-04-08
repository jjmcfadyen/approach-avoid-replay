function [startvals,fitvals,err] = cluster_fit_startrange(filename)
% fits model using a range of starting vlaues

tmp = load(filename);

model = tmp.model;
d = tmp.thisd;
numStart = tmp.numStart;
info = tmp.info;

[FILEPATH,NAME,EXT] = fileparts(filename);

% get function for this model
fun_model = @(x) sim_model(d,model,x);

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

if numStart^np > 100
    error(['Too many combinations of starting values to compute (' num2str(numStart^np) ')' ])
end

if size(startvals,1)==2
    startvals = combvec(startvals(1,:),startvals(2,:));
elseif size(startvals,1)==3
    startvals = combvec(startvals(1,:),startvals(2,:),startvals(3,:));
elseif size(startvals,1)==4
    startvals = combvec(startvals(1,:),startvals(2,:),startvals(3,:),startvals(4,:));
end

% optimise free parameters for each combination of starting values
optstart = nan(size(startvals,2),np+1);
for v = 1:size(startvals,2)
    if numStart>1
        disp(['--- Model ' model.name ': optimising parameters for start set ' num2str(v) ' of ' num2str(size(startvals,2)) '...'])
    end
    [x,err] = patternsearch(fun_model,startvals(:,v),[],[],[],[],lb,ub);
    optstart(v,:) = [err x'];
end

err = optstart(:,1)';
fitvals = optstart(:,2:end)';

% save
save(fullfile(FILEPATH,[NAME EXT]),'startvals','fitvals','err','d','model','info');

end