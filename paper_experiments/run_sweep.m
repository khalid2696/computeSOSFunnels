% P3: (rho_g, c) sweep with FIXED inlet/outlet regions -> Table 2 (and Fig. 4 curves).
% Requires ./results/default.mat from run_default.m.

clear inletScale runLog maxIter multiplierPolyDeg convergenceTolerance ...
      rhoInitialGuessConstant rhoInitialGuessExpCoeff captureSolverLog
load('./results/default.mat', 'fixedRegions');

drv_rho = [0.02 0.035 0.05 0.0625];   % rho_g  (outlet containment needs rho_g <= 0.0625)
drv_c   = [1 1.5 2 3];                % c      (inlet containment needs rho_g*e^c >= ~0.109)
drv_sweep = cell(numel(drv_rho), numel(drv_c));

for drv_i = 1:numel(drv_rho)
    for drv_j = 1:numel(drv_c)
        clear runLog
        rhoInitialGuessConstant = drv_rho(drv_i);
        rhoInitialGuessExpCoeff = drv_c(drv_j);
        maxIter = 4; convergenceTolerance = -inf; 
        multiplierPolyDeg = 6; drawFlag = 0;
        try
            run('step4_computeFunnel.m');
        catch drv_err
            fprintf('  (stopped: %s)\n', drv_err.message);
        end
        drv_sweep{drv_i, drv_j} = runLog;
        fprintf('rho_g=%.4g c=%.3g -> %s, iters=%d\n', drv_rho(drv_i), drv_c(drv_j), runLog.status, runLog.iters);
    end
end
save('./results/sweep.mat', 'drv_sweep', 'drv_rho', 'drv_c', 'fixedRegions');

fprintf('\nRatio floor sqrt(det M_T / det M_I) = %.4g\n', drv_sweep{1,1}.ratioFloor);

%% ---- LaTeX cells for Table 2: funnel vol | inlet-outlet ratio ----

load('./results/sweep.mat');
for i = 1:numel(drv_rho)
    row = sprintf('%.4g', drv_rho(i));
    for j = 1:numel(drv_c)
        L = drv_sweep{i,j};
        if strcmp(L.status,'infeasible_C') || drv_rho(i) >= 0.0625 - 1e-12, c = 'C';
        elseif L.iters == 0, c = '--';
        else, c = sprintf('%.3f / %.0f', L.funnelVol(end), L.inOutRatio); end
        row = [row ' & ' c];
    end
    fprintf('%s \\\\\n', row);
end