% Fig. 4: objective and funnel volume vs. alternation iteration (R1 data).
% Input : ./results/default.mat  (R1: maxIter = 20, convergenceTolerance = -Inf)
% Output: ./figures/alternation_funnel_chapter.png  + the numbers quoted in the text

load('./results/default.mat', 'runLog');

iterAxis = 0:runLog.iters;                       % iteration 0 = initial guess
J        = [runLog.J0, runLog.J];
vol      = [runLog.funnelVol0, runLog.funnelVol];
Jn       = J/J(1);                               % normalized by initial guess
voln     = vol/vol(1);

% ---- numbers for the text (check against the paragraph) ----
relDec   = -diff(runLog.J)./runLog.J(1:end-1);
iterStar = find(relDec < 0.05, 1) + 1;
[volMin, iMin] = min(runLog.funnelVol);
[volMax, iMax] = max(runLog.funnelVol(iterStar:end)); iMax = iMax + iterStar - 1;
fprintf('J strictly decreasing over all iterations: %d\n', all(diff(runLog.J) < 0));
fprintf('iter* = %d | captured = %.1f%% of 20-iteration decrease in J\n', ...
        iterStar, 100*(runLog.J0 - runLog.J(iterStar))/(runLog.J0 - runLog.J(end)));
fprintf('volume drop iter 1 -> iter*: %.1f%% | drop guess -> iter*: %.1f%%\n', ...
        100*(1 - runLog.funnelVol(iterStar)/runLog.funnelVol(1)), 100*(1 - runLog.funnelVol(iterStar)/runLog.funnelVol0));
fprintf('min volume at iter %d: vol(iter*) is %.1f%% above it\n', iMin, 100*(runLog.funnelVol(iterStar)/volMin - 1));
fprintf('largest volume after iter*: iter %d, %.1f%% above min\n', iMax, 100*(volMax/volMin - 1));

% ---- figure (single IEEE column: 3.5 in wide) ----
fig = figure('Units','inches','Position',[1 1 3.5 2.3],'Color','w');
hold on; box on; grid on;
plot(iterAxis, Jn,   '-o',  'LineWidth',1.2, 'MarkerSize',3, 'DisplayName','Objective $J_{\mathrm{iter}}/J_0$');
plot(iterAxis, voln, '--s', 'LineWidth',1.2, 'MarkerSize',3, 'DisplayName','Funnel volume (normalized)');
xline(iterStar, ':k', 'LineWidth',1.0, 'HandleVisibility','off');
text(iterStar + 0.3, 0.5, '$\mathrm{iter}^\star$', 'Interpreter','latex', 'FontSize',9, ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');
xlabel('Iteration, $\mathrm{iter}$', 'Interpreter','latex');
ylabel('Normalized value', 'Interpreter','latex');
xlim([0 runLog.iters]); ylim([0 1.05]);
legend('Interpreter','latex', 'Location','northeast', 'FontSize',8);
set(gca, 'FontSize',9, 'TickLabelInterpreter','latex');

if ~exist('./figures','dir'), mkdir('./figures'); end
exportgraphics(fig, './figures/alternation_funnel_chapter.png', 'Resolution', 600);
