clc; clearvars; close all

%% Add directories
addpath('../lib/');

%% Load the nominal trajectory and LQR gains
load('../precomputedData/nominalTrajectory.mat');
N = length(time_instances);

% -- Status message: Quantities at our disposal now -- %

% time_instances   : time horizon (sampled)        : 1 x N
% x_nom            : nominal state trajectory      : n_x x N
% u_nom            : nominal input tape            : n_u x N
% dynamicsFnHandle : the function handle of the system dynamics

%Specify the initial state perturbation
if ~exist('statePerturbation','var')
    statePerturbation = 0;
end

% Linear interpolation of inputs
u_func = @(t) interp1(time_instances, u_nom, t, 'linear', 'extrap');

%system dynamics
%dynamicsFnHandle is of the form @(x,u) cartpole_dynamics(x,u,cartPoleParameters);
cartpole_dynamics_odeFn = @(t,x) dynamicsFnHandle(x, u_func(t));

% Initial condition
x0 = x_nom(:,1) + statePerturbation;

% Time span
tspan = [time_instances(1), time_instances(end)];

% ODE solve
[t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, tspan, x0);
%[t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, time_instances, x0);
x_sol = x_sol'; t_sol = t_sol'; %to match the dimensions of x_nom

%% Plot state trajectories
figure;

subplot(2,2,1); hold on; grid on;
plot(t_sol, x_sol(1,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(1,:), 'k--', 'LineWidth', 1.75);
xlabel('t (s)'); ylabel('Cart Position (m)');

subplot(2,2,2); hold on; grid on;
plot(t_sol, x_sol(2,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(2,:), 'k--', 'LineWidth', 1.75);
xlabel('t (s)'); ylabel('Cart Velocity (m/s)');

subplot(2,2,3); hold on; grid on;
plot(t_sol, x_sol(3,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(3,:), 'k--', 'LineWidth', 1.75);
xlabel('t (s)'); ylabel('\theta (rad)');

subplot(2,2,4); hold on; grid on;
plot(t_sol, x_sol(4,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(4,:), 'k--', 'LineWidth', 1.75);
xlabel('t (s)'); ylabel('\theta dot (rad/s)');

sgtitle('Cart-Pole State Trajectories');

%% phase portraits
plotPhasePortraits(x_sol, x_nom, [1 3]);
drawnow;

%% input profiles
figure; grid on;
plot(time_instances, u_nom, 'LineWidth', 1.75);
xlabel('t (s)'); ylabel('F (N)');

%% Function definitions
function plotPhasePortraits(x_sol, x_nom, stateDims)
    % Plot all trajectories and nominal trajectory
    figure; hold on; grid on; axis equal;

    if nargin < 3
        stateDims = [1 3];
    end
    
    plot(x_nom(stateDims(1), :), x_nom(stateDims(2), :), 'k--', 'LineWidth', 2);
    plot(x_sol(stateDims(1), :), x_sol(stateDims(2), :), 'b-', 'LineWidth', 1.5);

    xlabel(['x_{', num2str(stateDims(1)), '}'])
    ylabel(['x_{', num2str(stateDims(2)), '}'])
    title('Monte Carlo Rollout Trajectories');
    legend('Nominal Trajectory','IVP solution','Location','best');
    %hold off;
end
