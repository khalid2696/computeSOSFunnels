% Fig. 2: representative funnel (R2 data = funnel output by Algorithm 1 at termination).
% Input : ./results/default_delta.mat (runLog from R2), ./precomputedData/nominalTrajectory.mat (x_nom)
% Output: ./figures/funnel_xy.png, ./figures/funnel_3d.png
close all

load('./results/representative.mat', 'runLog');
load('../precomputedData/nominalTrajectory.mat', 'x_nom');

P = runLog.P;  rho = runLog.rho;  M_I = runLog.M_I;  M_T = runLog.M_T;
N = size(x_nom, 2);
colI = [0.13 0.55 0.13];  colT = [0.80 0.15 0.15];  colF = [0.20 0.30 0.85];

projXY = @(M) inv([1 0; 0 1; 0 0]' * inv(M) * [1 0; 0 1; 0 0]);   % shadow of {x'Mx<=1} on (x,y)
th = linspace(0, 2*pi, 200);
ellXY = @(c, A) c + sqrtm(inv(A)) * [cos(th); sin(th)];            % boundary of {z'Az<=1} + c

% ---------------- (a) xy-projection ----------------
figA = figure('Units','inches','Position',[1 1 2.6 3.0],'Color','w'); hold on; box on; grid on; axis equal;
for k = 1:N
    e = ellXY(x_nom(1:2,k), projXY(P(:,:,k)/rho(k)));
    plot(e(1,:), e(2,:), '-', 'Color', [colF 0.6], 'LineWidth', 0.6, 'HandleVisibility','off');
end
eI = ellXY(x_nom(1:2,1),   projXY(M_I));
eT = ellXY(x_nom(1:2,end), projXY(M_T));
plot(x_nom(1,:), x_nom(2,:), '--k', 'LineWidth', 1.0, 'DisplayName', 'Nominal');
plot(eI(1,:), eI(2,:), '-.', 'Color', colI, 'LineWidth', 1.4, 'DisplayName', '$\mathcal{E}_I$');
plot(eT(1,:), eT(2,:), '-.', 'Color', colT, 'LineWidth', 1.4, 'DisplayName', '$\mathcal{E}_T$');
xlabel('$x$ [m]', 'Interpreter','latex'); ylabel('$y$ [m]', 'Interpreter','latex');
legend('Interpreter','latex', 'Location','southeast', 'FontSize',8);
set(gca, 'FontSize',9, 'TickLabelInterpreter','latex');

% ---------------- (b) 3D tube in (x, y, theta) ----------------
[Xs, Ys, Zs] = sphere(24);  U = [Xs(:) Ys(:) Zs(:)]';
ell3 = @(c, A) c + sqrtm(inv(A)) * U;                               % surface of {z'Az<=1} + c
figB = figure('Units','inches','Position',[1 1 3.4 3.0],'Color','w'); hold on; box on; grid on;
surf(reshape(SI(1,:), size(Xs)), reshape(SI(2,:), size(Xs)), reshape(SI(3,:), size(Xs)), ...
     'FaceColor', colI, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{E}_I$');

for k = 1:1:N
    S = ell3(x_nom(:,k), P(:,:,k)/rho(k));
    surf(reshape(S(1,:), size(Xs)), reshape(S(2,:), size(Xs)), reshape(S(3,:), size(Xs)), ...
         'FaceColor', colF, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility','off');
end
SI = ell3(x_nom(:,1), M_I);  ST = ell3(x_nom(:,end), M_T);
surf(reshape(ST(1,:), size(Xs)), reshape(ST(2,:), size(Xs)), reshape(ST(3,:), size(Xs)), ...
     'FaceColor', colT, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{E}_T$');
plot3(x_nom(1,:), x_nom(2,:), x_nom(3,:), '--k', 'LineWidth', 1.0, 'DisplayName', 'Nominal');
xlabel('$x$ [m]', 'Interpreter','latex'); ylabel('$y$ [m]', 'Interpreter','latex');
zlabel('$\theta$ [rad]', 'Interpreter','latex');
view(63, 25); camlight headlight; lighting gouraud;
legend('Interpreter','latex', 'Location','northeast', 'FontSize',8);
set(gca, 'FontSize',9, 'TickLabelInterpreter','latex');

if ~exist('./figures','dir'), mkdir('./figures'); end
exportgraphics(figA, './figures/funnel_xy.png', 'Resolution', 600);
exportgraphics(figB, './figures/funnel_3d.png', 'Resolution', 600);

% theta projections
load('../precomputedData/nominalTrajectory.mat', 'time_instances');
hw = arrayfun(@(k) sqrt(subsref(inv(P(:,:,k)/rho(k)), struct('type','()','subs',{{3,3}}))), 1:N);
hwI = sqrt(subsref(inv(M_I), struct('type','()','subs',{{3,3}})));
hwT = sqrt(subsref(inv(M_T), struct('type','()','subs',{{3,3}})));
t = time_instances;  th_nom = x_nom(3,:);

figB = figure('Units','inches','Position',[1 1 3.4 2.4],'Color','w'); hold on; box on; grid on;
errorbar(t(1),   th_nom(1),   hwI, 'Color',colI, 'LineWidth',1.4, 'CapSize',6, 'DisplayName','$\mathcal{E}_I$');
fill([t fliplr(t)], [th_nom+hw fliplr(th_nom-hw)], colF, 'FaceAlpha',0.2, 'EdgeColor','none', 'DisplayName','Funnel');
plot(t, th_nom, '--k', 'LineWidth',1.0, 'DisplayName','Nominal');
errorbar(t(end), th_nom(end), hwT, 'Color',colT, 'LineWidth',1.4, 'CapSize',6, 'DisplayName','$\mathcal{E}_T$');
xlabel('$t$ [s]', 'Interpreter','latex'); ylabel('$\theta$ [rad]', 'Interpreter','latex');
legend('Interpreter','latex', 'Location','best', 'FontSize',8);
set(gca, 'FontSize',9, 'TickLabelInterpreter','latex');
exportgraphics(figB, './figures/funnel_theta.png', 'Resolution', 600);