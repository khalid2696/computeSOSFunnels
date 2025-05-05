clc; clearvars; close all

%% Loading the time instances, nominal trajectory, nominal input and feedback control gain

load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/LQRGainsAndCostMatrices.mat');

% Load nominal trajectory and inputs
% x_nom: 3xN -> [x; y; theta]
% u_nom: 2xN -> [v; omega]
% K: 2x3xN -> Feedback gain matrices at each time step

%% Parameters
% Number of samples around each nominal point
num_samples = 50; 
% Define perturbation range around nominal trajectory
perturb_range = 0.5;  % Small perturbations in x, y, and theta

%% Start of analysis

% Number of time samples
N = length(time_instances);

% Extract trajectory and inputs
x0 = x_nom(1, :);
y0 = x_nom(2, :);
theta0 = x_nom(3, :);

% Generate perturbed states around nominal trajectory
[x_pert, y_pert, theta_pert] = deal([]);
[xdot_pert, ydot_pert, thetadot_pert] = deal([]);

for i = 1:N
    % Create small perturbations around the nominal state
    x_var = x0(i) + perturb_range * (2*rand(1,num_samples)-1);
    y_var = y0(i) + perturb_range * (2*rand(1,num_samples)-1);
    theta_var = theta0(i) + perturb_range * (2*rand(1,num_samples)-1);

    % Compute perturbed control using state feedback
    for j = 1:num_samples
        x_dev = [x_var(j); y_var(j); theta_var(j)] - x_nom(:,i); % Deviation
        u_pert = u_nom(:,i) - K(:,:,i) * x_dev; % Feedback control

        % Compute state derivatives
        xdot = u_pert(1) * cos(theta_var(j));
        ydot = u_pert(1) * sin(theta_var(j));
        thetadot = u_pert(2);

        % Store perturbed states and derivatives
        x_pert = [x_pert, x_var(j)];
        y_pert = [y_pert, y_var(j)];
        theta_pert = [theta_pert, theta_var(j)];

        xdot_pert = [xdot_pert, xdot];
        ydot_pert = [ydot_pert, ydot];
        thetadot_pert = [thetadot_pert, thetadot];
    end
end

% 3D Quiver plot showing flow field
figure; view(3);
quiver3(x_pert, y_pert, theta_pert, xdot_pert, ydot_pert, thetadot_pert, 'r', 'AutoScale', 'on', 'LineWidth', 1.2);
hold on;
plot3(x0, y0, theta0, 'b-', 'LineWidth', 2); % Nominal trajectory
xlabel('x'); ylabel('y'); zlabel('\theta');
title('State Flow Around Nominal Trajectory');
grid on; axis equal;
legend('Flow Field', 'Nominal Trajectory');

hold off;
