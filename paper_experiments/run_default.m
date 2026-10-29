% P2 + Fig. 4 baseline: default run. Defines fixedRegions used by ALL other experiments.
% Run this first, once.

clear fixedRegions inletScale runLog maxIter multiplierPolyDeg convergenceTolerance ...
      rhoInitialGuessConstant rhoInitialGuessExpCoeff captureSolverLog

rhoInitialGuessConstant = 0.05;   % rho_g
rhoInitialGuessExpCoeff = 1;      % c
maxIter                 = 20;     % j_max
multiplierPolyDeg       = 6;      % deg lambda_k
drawFlag                = 1;      % keep the plots for Fig. 2

run('./step4_computeFunnel.m');

%%
fixedRegions = struct('M_I', runLog.M_I, 'M_T', runLog.M_T);
if ~exist('./results','dir'), mkdir('./results'); end
save('./results/representative.mat', 'runLog', 'fixedRegions');

relDec = -diff(runLog.J)./runLog.J(1:end-1);
fprintf('\nDEFAULT: status=%s | iters=%d | converged=%d | regionsConvergent=%d\n', ...
        runLog.status, runLog.iters, runLog.converged, runLog.regionsConvergent);
fprintf('J history:         %s\n', mat2str(runLog.J, 5));
fprintf('rel. decrease:     %s\n', mat2str(relDec, 3));
fprintf('funnel vol history: %s\n', mat2str(runLog.funnelVol, 5));
