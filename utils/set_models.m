function models = set_models()

hightau = 3;
mintau = realmin;

models = struct();

m = 0;

% === FREE RANDOMNESS ============================================================================================================

% -----------------------------------------------
% Optimal two-path calculation (with free randomness)
% -----------------------------------------------

m = m+1;
models(m).name = 'optimal';

models(m).fun = 'twopath';
models(m).nThresh = 1;
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'threshold','tau'};
models(m).paraminfo.nFree = 2;

% -----------------------------------------------
% Only calculate path 1 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'path1';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path1'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'threshold','tau'};
models(m).paraminfo.nFree = 2;

% -----------------------------------------------
% Only calculate path 2 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'path2';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path2'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'threshold','tau'};
models(m).paraminfo.nFree = 2;

% -----------------------------------------------
% Only calculate one of the two paths (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'random';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'random'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'threshold','tau'};
models(m).paraminfo.nFree = 2;

% -----------------------------------------------
% Only calculate path 1 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'path1double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path1'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only calculate path 2 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'path2double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path2'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only calculate one of the two paths (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'randomdouble';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'random'; % which path is calculated
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only calculate rewarding path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'rewardinglearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'rewarding'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only calculate aversive path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'aversivelearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'aversive'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% === Q LEARNER ==============================================================================================

% -----------------------------------------------
% Two-path evaluation (with free randomness)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qoptimal';

models(m).fun = 'twopath';
models(m).nThresh = 1;
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 1 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qpath1';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path1'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 2 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qpath2';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path2'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate one of the two paths (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qrandom';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'random'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 1 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qpath1double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path1'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evaluate path 2 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qpath2double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path2'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evluate one of the two paths (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qrandomdouble';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'random'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evaluate rewarding path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qrewardinglearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'rewarding'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate aversive path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'Qaversivelearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'aversive'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% === Q LEARNER - NEGS ==============================================================================================

% -----------------------------------------------
% Two-path evaluation
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGoptimal';

models(m).fun = 'twopath';
models(m).nThresh = 1;
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 1 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGpath1';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path1'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 2 (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGpath2';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'path2'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate one of the two paths (one threshold)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGrandom';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'random'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate path 1 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGpath1double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path1'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evaluate path 2 (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGpath2double';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'path2'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evluate one of the two paths (two thresholds)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGrandomdouble';

models(m).fun = 'onepath';
models(m).nThresh = 2;
models(m).type = 'random'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.rewthreshold.val = NaN;
models(m).params.rewthreshold.start = 0;
models(m).params.rewthreshold.lb = -12;
models(m).params.rewthreshold.ub = 12;

models(m).params.lossthreshold.val = NaN;
models(m).params.lossthreshold.start = 0;
models(m).params.lossthreshold.lb = -12;
models(m).params.lossthreshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau','alpha'};
models(m).paraminfo.nFree = 4;

% -----------------------------------------------
% Only evaluate rewarding path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGrewardinglearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'rewarding'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% -----------------------------------------------
% Only evaluate aversive path (with learning)
% -----------------------------------------------

m = m+1;
models(m).name = 'QNEGaversivelearn';

models(m).fun = 'onepath';
models(m).nThresh = 1;
models(m).type = 'aversive'; % which path is calculated
models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
models(m).qlearner = true; % don't calculate value, just learn it from experience
models(m).neglearn = true; % store separate values for each position of odd rule per path

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = NaN;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = realmin;
models(m).params.tau.ub = Inf;

models(m).params.alpha.val = NaN;
models(m).params.alpha.start = 0.5;
models(m).params.alpha.lb = 0;
models(m).params.alpha.ub = 1;

models(m).paraminfo.names = {'threshold','tau','alpha'};
models(m).paraminfo.nFree = 3;

% === FIXED TAU (HIGHTAU) ==============================================================================================

% % -----------------------------------------------
% % Optimal two-path calculation (fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'optimalfixed';
% 
% models(m).fun = 'twopath';
% models(m).nThresh = 1;
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN; 
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau; % set to high value to make decision-making highly guided by value (minimal exploration)
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold'};
% models(m).paraminfo.nFree = 1;
% 
% % ------------------------------------------------
% % Only calculate path 1 (one threshold, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path1fixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'path1'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold'};
% models(m).paraminfo.nFree = 1;
% 
% % -----------------------------------------------
% % Only calculate path 2 (one threshold, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path2fixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'path2'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold'};
% models(m).paraminfo.nFree = 1;
% 
% % -----------------------------------------------
% % Only calculate one of the two paths (one threshold, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'randomfixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'random'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold'};
% models(m).paraminfo.nFree = 1;
% 
% % -----------------------------------------------
% % Only calculate path 1 (two thresholds, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path1doublefixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'path1'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate path 2 (two thresholds, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path2doublefixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'path2'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate one of the two paths (two thresholds, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'randomdoublefixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'random'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate rewarding path (with learning, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'rewardinglearnfixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'rewarding'; % which path is calculated
% models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).params.alpha.val = NaN;
% models(m).params.alpha.start = 0.5;
% models(m).params.alpha.lb = 0;
% models(m).params.alpha.ub = 1;
% 
% models(m).paraminfo.names = {'threshold','alpha'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate aversive path (with learning, fixed tau)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'aversivelearnfixed';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'aversive'; % which path is calculated
% models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = hightau;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = realmin;
% models(m).params.tau.ub = Inf;
% 
% models(m).params.alpha.val = NaN;
% models(m).params.alpha.start = 0.5;
% models(m).params.alpha.lb = 0;
% models(m).params.alpha.ub = 1;
% 
% models(m).paraminfo.names = {'threshold','alpha'};
% models(m).paraminfo.nFree = 2;
% 
% % === CONSTRAINED RANDOMNESS (MINTAU) ================================================================================
% 
% % -----------------------------------------------
% % Optimal two-path calculation (with constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'optimalconstrained';
% 
% models(m).fun = 'twopath';
% models(m).nThresh = 1;
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau; % make it highly likely to be value-guided (if EV of approaching is only 1 point higher than avoiding, there is a 95% chance of approaching)
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold','tau'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate path 1 (one threshold, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path1constrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'path1'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold','tau'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate path 2 (one threshold, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path2constrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'path2'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold','tau'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate one of the two paths (one threshold, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'randomconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'random'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'threshold','tau'};
% models(m).paraminfo.nFree = 2;
% 
% % -----------------------------------------------
% % Only calculate path 1 (two thresholds, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path1doubleconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'path1'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
% models(m).paraminfo.nFree = 3;
% 
% % -----------------------------------------------
% % Only calculate path 2 (two thresholds, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'path2doubleconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'path2'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
% models(m).paraminfo.nFree = 3;
% 
% % -----------------------------------------------
% % Only calculate one of the two paths (two thresholds, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'randomdoubleconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 2;
% models(m).type = 'random'; % which path is calculated
% models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.rewthreshold.val = NaN;
% models(m).params.rewthreshold.start = 0;
% models(m).params.rewthreshold.lb = -12;
% models(m).params.rewthreshold.ub = 12;
% 
% models(m).params.lossthreshold.val = NaN;
% models(m).params.lossthreshold.start = 0;
% models(m).params.lossthreshold.lb = -12;
% models(m).params.lossthreshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).paraminfo.names = {'rewthreshold','lossthreshold','tau'};
% models(m).paraminfo.nFree = 3;
% 
% % -----------------------------------------------
% % Only calculate rewarding path (with learning, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'rewardinglearnconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'rewarding'; % which path is calculated
% models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).params.alpha.val = NaN;
% models(m).params.alpha.start = 0.5;
% models(m).params.alpha.lb = 0;
% models(m).params.alpha.ub = 1;
% 
% models(m).paraminfo.names = {'threshold','tau','alpha'};
% models(m).paraminfo.nFree = 3;
% 
% % -----------------------------------------------
% % Only calculate aversive path (with learning, constrained randomness)
% % -----------------------------------------------
% 
% m = m+1;
% models(m).name = 'aversivelearnconstrained';
% 
% models(m).fun = 'onepath';
% models(m).nThresh = 1;
% models(m).type = 'aversive'; % which path is calculated
% models(m).learning = true; % false if they do the calculation, true if they learn the overall value of each path over time
% 
% % Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
% models(m).params.threshold.val = NaN;
% models(m).params.threshold.start = 0;
% models(m).params.threshold.lb = -12;
% models(m).params.threshold.ub = 12;
% 
% models(m).params.tau.val = NaN;
% models(m).params.tau.start = 0.1;
% models(m).params.tau.lb = mintau;
% models(m).params.tau.ub = Inf;
% 
% models(m).params.alpha.val = NaN;
% models(m).params.alpha.start = 0.5;
% models(m).params.alpha.lb = 0;
% models(m).params.alpha.ub = 1;
% 
% models(m).paraminfo.names = {'threshold','tau','alpha'};
% models(m).paraminfo.nFree = 3;

% === NULL =========================================================================================

% -----------------------------------------------
% Null model (preference for one option)
% -----------------------------------------------

m = m+1;
models(m).name = 'null';

models(m).fun = 'random';
models(m).nThresh = 1;
models(m).learning = false; % false if they do the calculation, true if they learn the overall value of each path over time

% Set fixed/free parameters (set to NaN to optimise it, otherwise sim_model.m uses defaults)
models(m).params.threshold.val = NaN;
models(m).params.threshold.start = 0;
models(m).params.threshold.lb = -12;
models(m).params.threshold.ub = 12;

models(m).params.tau.val = hightau;
models(m).params.tau.start = 0.1;
models(m).params.tau.lb = mintau;
models(m).params.tau.ub = Inf;

models(m).paraminfo.names = {'threshold'};
models(m).paraminfo.nFree = 1;

end