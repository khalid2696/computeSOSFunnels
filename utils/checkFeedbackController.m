clc; close all 
clearvars;

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

cartpole_dynamics = dynamicsFnHandle;

%% Feedback controller synthesis (need to be run only during design)

% nx = size(x_nom, 1); nu = size(u_nom, 1); N = length(time_instances);
% x = createSymbolicVector('x', nx); %state vector
% u = createSymbolicVector('u', nu); %input vector
% % [x1, x2, x3, x4] -- [px; vx; theta; omega]
% % u1 -- F
% 
% f = dynamicsFnHandle(x, u);
% 
% Q = 10*eye(nx);
% R = 0.1;
% 
% [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom(:,end), u_nom(:,end), Q, R);
% 
% keyboard
%% Forward rollout

startTimeOffset = N-1;
x_nom = x_nom(:,end-startTimeOffset:end);
u_nom = u_nom(:,end-startTimeOffset:end);
time_instances = time_instances(end-startTimeOffset:end);

% Initial condition
statePerturbation = 0; %5e-3
%x0 = x_nom(:,1) + statePerturbation*ones(size(x_nom,1),1);
x0 = x_nom(:,1) + statePerturbation*rand(size(x_nom,1),1);
xf = x_nom(:,end);

% Linear interpolation of inputs
u_func = @(t) interp1(time_instances, u_nom, t, 'spline', 'extrap');

%system dynamics
%dynamicsFnHandle is of the form @(x,u) cartpole_dynamics(x,u,cartPoleParameters);
cartpole_dynamics_odeFn = @(t,x) cartpole_dynamics(x, u_func(t));

% Time span
tspan = [time_instances(1), time_instances(end)];

% ODE solve
[t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, tspan, x0);
%[t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, time_instances, x0);
x_sol = x_sol'; t_sol = t_sol'; %to match the dimensions of x_nom

% %using the same scheme as how the traj optimisation was solved (RK4)
% t_sol = time_instances; dt = time_instances(2) - time_instances(1);
% x_sol = NaN(size(x_nom));
% 
% x_sol(:,1) = x0;
% for k = 1:length(time_instances)-1
% 
%     xk = x_sol(:,k);
%     uk = u_nom(k);
% 
%     % %trapezoidal
%     % f_k = cartpole_dynamics(xk, uk);
%     % f_k_next = cartpole_dynamics(xk + f_k*dt, uk);
%     % x_next = xk + (dt/2) * (f_k + f_k_next);
% 
%     %RK4 integration    
%     k1 = cartpole_dynamics(xk, uk);
%     k2 = cartpole_dynamics(xk + 0.5 * dt * k1, uk);
%     k3 = cartpole_dynamics(xk + 0.5 * dt * k2, uk);
%     k4 = cartpole_dynamics(xk + dt * k3, uk);
%     x_next = xk + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4); % Update state
% 
%     x_sol(:,k+1) = x_next;
% end

%% LQR based stabilisation phase

force_limits = [-12 12];

%LQR gain for stabilising at top
%K = [-3.1623   -5.6172   44.2219   10.2343];
K = [-10.0000  -16.2819   91.7720   22.6933];

%LQR gain for stabilising at hanging
% K = 10*[0.0707    0.3946   -0.0058   -0.0057];

stabilizeTimeHorizon = 7;
samplingRate = 100;
numTimeSteps = samplingRate*stabilizeTimeHorizon;
dt = 1/samplingRate;
T = linspace(time_instances(end),time_instances(end)+stabilizeTimeHorizon,numTimeSteps);

x_stabilize = NaN(4,numTimeSteps);
u_stabilize = NaN(1,numTimeSteps);
x_stabilize(:,1) = x_sol(:,end); %getting from previous forward integration

for k = 1:length(T)-1

    xk = x_stabilize(:,k);
    uk = K*(xf - xk);
    
    %force clipping
    uk = min(max(uk, force_limits(1)), force_limits(2));

    k1 = cartpole_dynamics(xk, uk);
    k2 = cartpole_dynamics(xk + 0.5 * dt * k1, uk);
    k3 = cartpole_dynamics(xk + 0.5 * dt * k2, uk);
    k4 = cartpole_dynamics(xk + dt * k3, uk);

    x_next = xk + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4); % Update state

    x_stabilize(:,k+1) = x_next;
    u_stabilize(:,k) = uk;
end

u_stabilize(:,end) = uk; %padding to maintain same array length
%% Plot state trajectories
figure;

subplot(2,2,1); hold on; grid on;
plot(t_sol, x_sol(1,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(1,:), 'k--', 'LineWidth', 1.75);
plot(T,x_stabilize(1,:), 'g:','LineWidth', 1.75);
xlabel('t (s)'); ylabel('Cart Position (m)');

subplot(2,2,2); hold on; grid on;
plot(t_sol, x_sol(2,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(2,:), 'k--', 'LineWidth', 1.75);
plot(T,x_stabilize(2,:), 'g:','LineWidth', 1.75);
xlabel('t (s)'); ylabel('Cart Velocity (m/s)');

subplot(2,2,3); hold on; grid on;
plot(t_sol, x_sol(3,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(3,:), 'k--', 'LineWidth', 1.75);
plot(T,x_stabilize(3,:), 'g:','LineWidth', 1.75);
xlabel('t (s)'); ylabel('\theta (rad)');

subplot(2,2,4); hold on; grid on;
plot(t_sol, x_sol(4,:), 'LineWidth', 1.5);
plot(time_instances, x_nom(4,:), 'k--', 'LineWidth', 1.75);
plot(T,x_stabilize(4,:), 'g:','LineWidth', 1.75);
xlabel('t (s)'); ylabel('\theta dot (rad/s)');

sgtitle('Cart-Pole State Trajectories');

%% phase portraits
plotPhasePortraits(x_sol, x_nom, x_stabilize, [1 3]);
drawnow;

%% input profiles
figure; grid on; hold on
plot(time_instances, u_nom, 'LineWidth', 1.75);
plot(T,u_stabilize, 'g:','LineWidth', 1.75);
xlabel('t (s)'); ylabel('F (N)');
title('Input profile');
drawnow
    
%% Function definitions

function sym_vector = createSymbolicVector(vecSymbol, n)
    % Generate variable names dynamically
    varNames = arrayfun(@(i) [vecSymbol, num2str(i)], 1:n, 'UniformOutput', false);
    
    % Define these as symbolic variables
    syms(varNames{:});
    
    % Construct a symbolic vector from these variables
    sym_vector = sym(varNames, 'real');

    sym_vector = sym_vector'; %column matrix by convention -- nx1
end

%an infinite horizon LQR to stabilize the closed loop system at the terminal state
function [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom_f, u_nom_f, Q, R)

    A_symb = jacobian(f, x);
    B_symb = jacobian(f, u);
    
    %Get the linearised matrices (Jacobians) at terminal state
    A_f = double(subs(A_symb, [x; u], [x_nom_f; u_nom_f]));
    B_f = double(subs(B_symb, [x; u], [x_nom_f; u_nom_f]));
    
    [K_f,S_f,CLP] = lqr(A_f,B_f,Q,R)
end


function plotPhasePortraits(x_sol, x_nom, x_hold, stateDims)
    % Plot all trajectories and nominal trajectory
    figure; hold on; grid on; axis equal;

    if nargin < 4
        stateDims = [1 3];
    end
    
    plot(x_nom(stateDims(1), :), x_nom(stateDims(2), :), 'k--', 'LineWidth', 2);
    plot(x_sol(stateDims(1), :), x_sol(stateDims(2), :), 'b-', 'LineWidth', 1.5);

    if nargin > 2
        plot(x_hold(stateDims(1), :), x_hold(stateDims(2), :), 'g:', 'LineWidth', 1.5);
    end

    xlabel(['x_{', num2str(stateDims(1)), '}'])
    ylabel(['x_{', num2str(stateDims(2)), '}'])
    title('IVP Rollout Trajectory');
    legend('Nominal Trajectory','IVP solution','LQR stabilisation', 'Location','best');
    %hold off;
end
