% P5: per-iteration timing and problem size -> Table 3.
% Copy into the cart-pole and quadrotor folders and set drv_deg accordingly:
%   unicycle: [2 4 6] | cart-pole: [2 4 6] | quadrotor: [2]
% Requires ./results/default.mat (run run_default.m in each folder first).

clear inletScale runLog maxIter multiplierPolyDeg convergenceTolerance ...
      rhoInitialGuessConstant rhoInitialGuessExpCoeff captureSolverLog
load('./results/default.mat', 'fixedRegions');

drv_deg  = [2];
drv_name = 'quadrotor';
rhoInitialGuessConstant = 1e-3;   % rho_g
rhoInitialGuessExpCoeff = 3;      % c
drv_tim  = struct('deg', {}, 'runLog', {});

% warm-up run (license checkout, JIT) -- discarded
maxIter = 1; multiplierPolyDeg = drv_deg(1); drawFlag = 0;
try run('step4_computeFunnel.m'); catch, end

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
%%
clear convergenceTolerance   % restore default (5e-2) for later runs
save(sprintf('./results/timing_%s.mat', drv_name), 'drv_tim');

%%
% Table 3: prints LaTeX rows from the timing results of all systems,
% checks measured PSD sizes against the analytic count, and prints the
% MOSEK share of wall time (for the "1--4%" sentence).
% Edit the paths below to where each system's results/timing_<name>.mat lives.

file = {'Quadrotor',  './results/timing_quadrotor.mat';};

rows = {};  share = [];
for f = 1:size(file,1)
    if ~isfile(file{f,2}), fprintf('%% missing: %s\n', file{f,2}); continue; end
    S = load(file{f,2}, 'drv_tim');
    for i = 1:numel(S.drv_tim)
        L   = S.drv_tim(i).runLog;  deg = S.drv_tim(i).deg;
        n   = L.n;  K = L.K;
        d   = max(2, (deg + 2)/2);                 % half-degree of the invariance constraint
        g   = nchoosek(n + d, d);                  % Gram size of each invariance block
        if L.iters < 1
            rows{end+1} = sprintf('%s & %d & %d & %d & %d & \\multicolumn{3}{c}{initial guess infeasible} \\\\', ...
                                  file{f,1}, n, K, deg, g);
            continue
        end
        % analytic scalarized PSD entries (dense Gram basis for invariance blocks)
        predL = K*g*(g+1)/2 + K*n*(n+1)/2       + 2*(n+1)*(n+2)/2 + 2;
        predV = K*g*(g+1)/2 + (K+1)*n*(n+1)/2   + (K+1) + 2*(n+1)*(n+2)/2 + 1;
        measL = L.sizeL(1,4);  measV = L.sizeV(1,4);
        if measL ~= predL || measV ~= predV
            fprintf('%% NOTE %s deg=%d: measured PSD entries %d/%d vs analytic %d/%d (SOSTOOLS pruned monomials?)\n', ...
                    file{f,1}, deg, measL, measV, predL, predV);
        end
        tL = mean(L.tL); tV = mean(L.tV); cL = mean(L.cpuL); cV = mean(L.cpuV);
        share = [share, cL/tL, cV/tV];
        rows{end+1} = sprintf('%s & %d & %d & %d & %d & %d / %d & %.3g (%.2g) & %.3g (%.2g) \\\\', ...
                              file{f,1}, n, K, deg, g, measL, measV, tL, cL, tV, cV); %#ok<SAGROW>
    end
end

fprintf('\n%% ---- paste into Table 3 ----\n');
fprintf('%s\n', rows{:});
fprintf('\n%% MOSEK share of wall time: %.1f%% to %.1f%%\n', 100*min(share), 100*max(share));