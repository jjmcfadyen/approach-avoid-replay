function hybrid = makeHybrid(MF,MB)
% function hybrid = makeHybrid(MF,MB)
% MF is the model-free model structure
% MB is the model-based model structure

hybrid = [];
hybrid.plantype = 'hybrid';
hybrid.MB = MB;
hybrid.MF = MF;

% Add replay parameters
hybrid.params.replayalpha.val = NaN;
hybrid.params.replayalpha.start = 0.9;
hybrid.params.replayalpha.lb = 0;
hybrid.params.replayalpha.ub = 1;

hybrid.params.gain.val = NaN;
hybrid.params.gain.start = 0.001;
hybrid.params.gain.lb = 0;
hybrid.params.gain.ub = 12;

% Unpack other model parameters
for m = {'MB','MF'}
    if ~isempty(hybrid.(m{1}))
        pnames = fieldnames(hybrid.(m{1}).params);
        for p = 1:length(pnames)
            if strcmp(pnames{p},'tau') || strcmp(pnames{p},'alpha')
                if ~isfield(hybrid.params,'tau') ||  ~isfield(hybrid.params,'alpha')
                    hybrid.params.(erase(pnames{p},m{1})) = hybrid.(m{1}).params.(pnames{p}); % only need one inverse temperature / learning rate parameter
                end
            else
                hybrid.params.([m{1} '_' pnames{p}]) = hybrid.(m{1}).params.(pnames{p});
            end
        end
    end
end
if any(ismember(fieldnames(hybrid.params),'MB_threshold')) && any(ismember(fieldnames(hybrid.params),'MF_threshold')) % only need threshold for model BASED planning
    hybrid.params.threshold = hybrid.params.MB_threshold;
    hybrid.params = rmfield(hybrid.params,'MB_threshold');
    hybrid.params = rmfield(hybrid.params,'MF_threshold');
end

if ~isempty(MB) && ~isempty(MF)
    if strcmp(MF.name,MB.name) % only need ONE forgetting parameter if models are the same
        if any(ismember(fieldnames(hybrid.params),'MB_forgetting')) && any(ismember(fieldnames(hybrid.params),'MF_forgetting')) % only need threshold for model BASED planning
            hybrid.params.forgetting = hybrid.params.MF_forgetting;
            hybrid.params = rmfield(hybrid.params,'MB_forgetting');
            hybrid.params = rmfield(hybrid.params,'MF_forgetting');
        end
    end
end

end