% P4: validity of the certificate under the TRUE dynamics, for inlet scales s -> Fig. 3(a).
% Requires ./results/default.mat from run_default.m.

clear inletScale runLog maxIter multiplierPolyDeg convergenceTolerance ...
      rhoInitialGuessConstant rhoInitialGuessExpCoeff captureSolverLog
load('./results/default.mat', 'fixedRegions');

drv_s   = [0.5 1 1.5];     % inlet scales: E_I -> s*E_I
drv_Ms  = 1000;            % boundary samples per knot (true dynamics)
drv_MsP = 50;              % boundary samples per knot (Taylor polynomial, sanity check; slow)
drv_val = struct('s', {}, 'runLog', {}, 'viol', {}, 'maxExcess', {}, 'violPoly', {});

for drv_i = 1:numel(drv_s)
    clear runLog
    inletScale              = drv_s(drv_i);
    rhoInitialGuessConstant = 0.05;
    rhoInitialGuessExpCoeff = 1 + 2*max(0, log(inletScale));  % enlarge guess inlet by s^2 only when s > 1
    maxIter = 20; multiplierPolyDeg = 6; drawFlag = 0;
    try
        run('step4_computeFunnel.m');
    catch drv_err
        fprintf('  (stopped: %s)\n', drv_err.message);
    end
    drv_val(drv_i).s      = inletScale;
    drv_val(drv_i).runLog = runLog;
    if runLog.iters >= 1
        [drv_val(drv_i).viol, drv_val(drv_i).maxExcess, drv_val(drv_i).violPoly] = ...
            checkBoundaryValidity(runLog, x_nom, u_nom, K, time_instances, deviationDynamics, xbar, ...
                                  @unicycleDynamics, drv_Ms, drv_MsP);
        fprintf('s=%.2f: mean violation (true) = %.4f | max over knots = %.4f | Taylor sanity max = %.4f\n', ...
                inletScale, mean(drv_val(drv_i).viol, 'omitnan'), max(drv_val(drv_i).viol), max(drv_val(drv_i).violPoly));
    else
        fprintf('s=%.2f: no certificate (%s)\n', inletScale, runLog.status);
    end
end
clear inletScale   % IMPORTANT: do not leak into later runs
save('./results/validity.mat', 'drv_val');
