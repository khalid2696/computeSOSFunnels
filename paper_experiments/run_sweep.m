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
        maxIter = 20; multiplierPolyDeg = 6; drawFlag = 0;
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

% ---- LaTeX cells for Table 2: funnel vol / inlet-outlet ratio / iterations (* = hit j_max) ----
fprintf('\nRatio floor sqrt(det M_T / det M_I) = %.4g\n', drv_sweep{1,1}.ratioFloor);
for drv_i = 1:numel(drv_rho)
    drv_row = sprintf('%.4g', drv_rho(drv_i));
    for drv_j = 1:numel(drv_c)
        L = drv_sweep{drv_i, drv_j};
        if strcmp(L.status, 'infeasible_C')
            drv_cell = 'C';
        elseif strcmp(L.status, 'infeasible_I') || L.iters == 0
            drv_cell = '--';
        else
            drv_star = ''; if ~L.converged, drv_star = '$^*$'; end
            drv_cell = sprintf('%.3g / %.3g / %d%s', L.funnelVol(end), L.inOutRatio, L.iters, drv_star);
        end
        drv_row = [drv_row ' & ' drv_cell]; %#ok<AGROW>
    end
    fprintf('%s \\\\\n', drv_row);
end
