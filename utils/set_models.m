function [models,hybrids] = set_models()

hightau = 3;
mintau = realmin;

models = struct();

m = 0;
for plantype = {'MB','MF','bias'} % model-based = using known probabilities, model-free = mapping value to actions, bias = null model with preference to approach/avoid

    switch plantype{1}
        case 'MB'
            CALC = {'calculation','experience'};
        case 'MF'
            CALC = {'actions'}; % only actions are considered by MF agent
        case 'bias'
            CALC = {'NA'};
    end

    for calctype = CALC

        switch calctype{1}
            case 'experience'
                QTYPE = {'perpath','perneg'}; % how many values to cache (one per path, or 3 per path depending on odd rule position)
            case 'calculation'
                QTYPE = {'NA'};
            case 'actions'
                QTYPE = {'peraction','perprob','perneg','pernegprob'}; % how many states to map actions to
            otherwise
                QTYPE = {'NA'};
        end
            
        for qtype = QTYPE
            
            switch plantype{1}
                case 'MB'
                    PATHTYPE = {'twopath','path1','path2','random','rewarding','aversive'}; % which paths are taken into consideration
                case 'MF'
                    PATHTYPE = {'NA'}; % only actions are considered by MF agent
                otherwise
                    PATHTYPE = {'NA'};
            end

            for pathtype = PATHTYPE

                switch pathtype{1}
                    case {'path1','path2','random'}
                        THRESHTYPE = {'single','double'};
                    case {'twopath','rewarding','aversive'}
                        THRESHTYPE = {'single'};
                    otherwise
                        THRESHTYPE = {'NA'};
                end
                if strcmp(plantype{1},'bias')
                    THRESHTYPE = {'single'};
                end

                for threshtype = THRESHTYPE

                    for forgetting = false %[true false]

%                         if forgetting
%                             forgetname = ['withforgetting'];
%                         else
%                             forgetname = ['noforgetting'];
%                         end
                        forgetname = [];

                        modelname = [plantype{1} '_' calctype{1} '_' qtype{1} '_' pathtype{1} '_' threshtype{1} '_' forgetname];
                        modelname = erase(modelname,'_NA');
    
                        m = m+1;
                        models(m).name = modelname;
    
                        % model family
                        models(m).plantype   = plantype{1};     % MB or MF
                        models(m).calctype   = calctype{1};     % mental calculation or learn from experience
                        models(m).qtype      = qtype{1};        % if learn from experience, then how many values to cache
                        models(m).pathtype   = pathtype{1};     % which path(s) to appraise
                        models(m).threshtype = threshtype{1};   % if single path, whether to use one or two thresholds, depending on sign (pos or neg) of path
    
                        % threshold parameter
                        if strcmp(models(m).plantype,'MB') || strcmp(models(m).plantype,'bias')
                            switch models(m).threshtype
                                case 'single'
                                    models(m).params.threshold.val = NaN;
                                    models(m).params.threshold.start = 0;
                                    models(m).params.threshold.lb = -12;
                                    models(m).params.threshold.ub = 12;
                                case 'double'
                                    models(m).params.rewthreshold.val = NaN;
                                    models(m).params.rewthreshold.start = 0;
                                    models(m).params.rewthreshold.lb = -12;
                                    models(m).params.rewthreshold.ub = 12;
        
                                    models(m).params.lossthreshold.val = NaN;
                                    models(m).params.lossthreshold.start = 0;
                                    models(m).params.lossthreshold.lb = -12;
                                    models(m).params.lossthreshold.ub = 12;
                            end
                        end
    
                        % learning rate parameter
                        if strcmp(models(m).plantype,'MF') || strcmp(models(m).calctype,'experience') || any(ismember(models(m).pathtype,{'rewarding','aversive'}))
                            models(m).params.alpha.val = NaN;
                            models(m).params.alpha.start = 0.5;
                            models(m).params.alpha.lb = 0;
                            models(m).params.alpha.ub = 1;
                        end
    
                        % inverse temperature parameter
                        models(m).params.tau.val = NaN;
                        models(m).params.tau.start = 0.1;
                        models(m).params.tau.lb = realmin;
                        models(m).params.tau.ub = Inf;

                        % forgetting parameter
                        if forgetting && ~strcmp(models(m).qtype,'NA')
                            models(m).params.forgetting.val = NaN;
                        else
                            models(m).params.forgetting.val = 0; % no forgetting (needs to be 0 for any calculation models)
                        end
                        models(m).params.forgetting.start = 0.001;
                        models(m).params.forgetting.lb = 0-realmin;
                        models(m).params.forgetting.ub = 1+realmin;
                    end
                end
            end
        end
    end
end

%% Hybrid models

% Pick apart different path selection strategies
mnames = extractfield(models,'name');

pathselection = {'twopath','path1_single','path2_single','random_single','path1_double','path2_double','random_double','rewarding','aversive'};
pathappraisal = {'calculation','experience_perpath','experience_perneg'};
pathlearning = {'experience_perpath','experience_perneg'};

% Create hybrid models, where a certain path appraisal strategy is used alongside a habitual learning of path values
cc = 0;
hybrids = [];
for strategy = pathselection
    for appraisal = pathappraisal % how does the model-based system compute the value of each path?
        for learning = pathlearning % how does the model-free system learn the value of each path?
            for teaching = {'mb-teach-mf','mf-teach-mb'} % which system does the updating?

                MB = contains(mnames,strategy{1}) & contains(mnames,appraisal{1});
                MF = contains(mnames,strategy{1}) & contains(mnames,learning{1});
        
                tmp = makeHybrid(models(MF),models(MB));
                tmp.name = ['Hybrid_' appraisal{1} '_' strategy{1} '_' learning{1}];
                tmp.gaintype = 'new-old';
                tmp.teachtype = teaching{1};
    
                cc = cc+1;
                if isempty(hybrids)
                    hybrids = tmp;
                else
                    hybrids(cc) = tmp;
                end
            end
        end
    end
end

% Make sure the model-based system has no forgetting
for m = 1:length(hybrids)

    hybrids(m).params.MB_forgetting.val = 0;
    hybrids(m).MB.params.forgetting.val = 0;

end

end