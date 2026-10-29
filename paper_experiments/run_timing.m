% P5: per-iteration timing and problem size -> Table 3.
% Copy into the cart-pole and quadrotor folders and set drv_deg accordingly:
%   unicycle: [2 4 6] | cart-pole: [2 6] | quadrotor: [2]
% Requires ./results/default.mat (run run_default.m in each folder first).

clear inletScale runLog maxIter multiplierPolyDeg convergenceTolerance ...
      rhoInitialGuessConstant rhoInitialGuessExpCoeff captureSolverLog
load('./results/default.mat', 'fixedRegions');

drv_deg  = [2 4 6];
drv_name = 'unicycle';
drv_tim  = struct('deg', {}, 'runLog', {});

% warm-up run (license checkout, JIT) -- discarded
maxIter = 1; multiplierPolyDeg = drv_deg(end); drawFlag = 0;
try, run('step4_computeFunnel.m'); catch, end

for drv_i = 1:numel(drv_deg)
    clear runLog
    multiplierPolyDeg    = drv_deg(drv_i);
    maxIter              = 3;
    convergenceTolerance = -Inf;      % force exactly 3 iterations for timing
    drawFlag             = 0;
    try
        run('step4_computeFunnel.m');
    catch drv_err
        fprintf('  (stopped: %s)\n', drv_err.message);
    end
    drv_tim(drv_i).deg    = drv_deg(drv_i);
    drv_tim(drv_i).runLog = runLog;
    if runLog.iters >= 1
        fprintf(['%s deg=%d | K=%d | L: wall %.3f s, MOSEK %.3f s | V: wall %.3f s, MOSEK %.3f s | ' ...
                 'L size [cons scal psd scalarized] = %s | V size = %s\n'], ...
                drv_name, drv_deg(drv_i), runLog.K, mean(runLog.tL), mean(runLog.cpuL), ...
                mean(runLog.tV), mean(runLog.cpuV), mat2str(runLog.sizeL(1,:)), mat2str(runLog.sizeV(1,:)));
    else
        fprintf('%s deg=%d: %s\n', drv_name, drv_deg(drv_i), runLog.status);
    end
end
clear convergenceTolerance   % restore default (1e-2) for later runs
save(sprintf('./results/timing_%s.mat', drv_name), 'drv_tim');
