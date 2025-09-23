clc; close all 
clearvars;

%% Add directories
addpath('../lib/');

%% Load the nominal trajectory and LQR gains
load('../precomputedData/nominalTrajectory.mat');
%N = length(time_instances);

% -- Status message: Quantities at our disposal now -- %

% time_instances   : time horizon (sampled)        : 1 x N
% x_nom            : nominal state trajectory      : n_x x N
% u_nom            : nominal input tape            : n_u x N
% dynamicsFnHandle : the function handle of the system dynamics

%Specify the initial state perturbation
if ~exist('intialStatePerturbation','var')
    intialStatePerturbation = 0;
end

numSamples = 500;

%% Initialisation parameters
nx = 4; % state-dimensionality 
nu = 1; % input-dimensionality

% Initial condition
intialStatePerturbation = 1e-3; %1e-2
x0 = -intialStatePerturbation/2*ones(nx,1) + intialStatePerturbation*rand(nx,1);
%x0 = zeros(4,1);
xf = [0 0 pi 0]';

%% Feedback controller synthesis (need to be run only during design)

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
% [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, xf, 0, Q, R);
% 
% keyboard

%% Control strategy parameters

%desired final state
controlStrategyParameters.xf = xf;
%gains
controlStrategyParameters.k_s = 0.75;
controlStrategyParameters.K = [-10.0000  -16.2819   91.7720   22.6933];
%switching control strategy
controlStrategyParameters.modeSwitchingAngle = 0.8*pi;
controlStrategyParameters.u0 = 5; %some initial input to warmstart the swing-up
controlStrategyParameters.initialImpulseTime = 0.1; %time duration of warmstart
%control saturations
controlStrategyParameters.force_limits = [-12 12];

%% Forward rollout

%system dynamics
cartpole_dynamics_odeFn = @(t,x) cartpole_dynamics(t, x, cartPoleParameters, controlStrategyParameters);

% Time span
N = 100; %number of time instances (discretisation)
%tspan = [0 15];
tspan = linspace(0, 10, N);

trajectories = cell(numSamples, 1);
input_profiles = cell(numSamples, 1);

for i = 1:numSamples
    x0 = -intialStatePerturbation/2*ones(nx,1) + intialStatePerturbation*rand(nx,1);
    
    % ODE solve
    [t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, tspan, x0);
    %[t_sol, x_sol] = ode45(cartpole_dynamics_odeFn, time_instances, x0);
    x_sol = x_sol'; t_sol = t_sol'; %to match the dimensions of x_nom
    
    % get the input profile after solving for the nominal trajectory
    u_sol = NaN(1,length(t_sol));
    for k = 1:length(t_sol)
        if t_sol(k) < controlStrategyParameters.initialImpulseTime
            u_sol(k) = controlStrategyParameters.u0;
        else
            u_sol(k) = swingUpEnergyControlLaw(x_sol(:,k), controlStrategyParameters);
        end
    end

    trajectories{i} = x_sol;
    input_profiles{i} = u_sol;
end

%% Plot state trajectories
figure;

subplot(2,2,1); hold on; grid on;
for i = 1:numSamples
    x_sol = trajectories{i};
    plot(t_sol, x_sol(1,:), 'b', 'LineWidth', 1.5);
end
xlabel('t (s)'); ylabel('Cart Position (m)');

subplot(2,2,2); hold on; grid on;
for i = 1:numSamples
    x_sol = trajectories{i};
    plot(t_sol, x_sol(2,:), 'b', 'LineWidth', 1.5);
end
xlabel('t (s)'); ylabel('Cart Velocity (m/s)');

subplot(2,2,3); hold on; grid on;
for i = 1:numSamples
    x_sol = trajectories{i};
    plot(t_sol, rad2deg(x_sol(3,:)), 'b', 'LineWidth', 1.5);
end
yline(rad2deg(controlStrategyParameters.modeSwitchingAngle), 'm--', 'Mode-switch');
yline(180, 'k--', 'Upright');
xlabel('t (s)'); ylabel('\theta (rad)');

subplot(2,2,4); hold on; grid on;
for i = 1:numSamples
    x_sol = trajectories{i};
    plot(t_sol, x_sol(4,:), 'b', 'LineWidth', 1.5);
end
xlabel('t (s)'); ylabel('\theta dot (rad/s)');

sgtitle('Cart-Pole State Trajectories');

%% phase portraits
plotPhasePortraits(trajectories, intialStatePerturbation, [1 3]);

%% input profiles
figure; grid on; hold on
for i = 1:numSamples
    u_sol = input_profiles{i};
    plot(t_sol, u_sol, 'k-.', 'LineWidth', 1.75);
end
xlabel('t (s)'); ylabel('F (N)');
title('Input profile');
drawnow

return

%% animation
%animate_cartpole_trajectory(t_sol, x_sol, cartPoleParameters, 1); %last argument is saveVideoFlag

%% save the nominal trajectory and feedforward input
time_instances = t_sol;
x_nom = x_sol;
u_nom = u_sol;

save('../precomputedData/nominalTrajectory.mat', 'time_instances', 'x_nom', 'u_nom', 'nom_trajCost', 'dynamicsFnHandle', 'cartPoleParameters');

disp('Saved the nominal trajectory and nominal inputs to a file!');
disp(' ');
%% Function definitions

function u = swingUpEnergyControlLaw(x, controlStrategyParameters)
    
    force_limits = controlStrategyParameters.force_limits;
    k_s = controlStrategyParameters.k_s;
    K = controlStrategyParameters.K;
    modeSwitchingAngle = controlStrategyParameters.modeSwitchingAngle;

    theta = x(3);
    omega = x(4);

    %final state
    xf = controlStrategyParameters.xf;
    
    %swing up control law
    if abs(theta) < modeSwitchingAngle
        u = -k_s*omega*cos(theta);
    else
        u = K*(xf - x);
    end

    %force clipping
    u = min(max(u, force_limits(1)), force_limits(2));

end

function f = cartpole_dynamics(t, x, cartPoleParameters, controlStrategyParameters)
    % Numerical evaluation of cartpole dynamics 
    
    %extract parameters
    M = cartPoleParameters.M; m = cartPoleParameters.m;
    L = cartPoleParameters.L; g = cartPoleParameters.g;

    % Extract states
    p_x = x(1);
    v_x = x(2);
    theta = x(3);
    omega = x(4);
    
    %swingup control law
    if t < controlStrategyParameters.initialImpulseTime %for some initial time
        F = controlStrategyParameters.u0; %some initial force to warm-start the swing up
    else
        F = swingUpEnergyControlLaw(x, controlStrategyParameters);
    end

    % Define trigonometric functions
    s_theta = sin(theta);
    c_theta = cos(theta);
    
    % Common denominator
    %denom = M + m*(1 - c_theta^2);
    denom = M + m*s_theta^2;
    
    % State derivatives
    f = [
        v_x;
        (F + m*L*omega^2*s_theta + m*g*s_theta*c_theta) / denom;
        omega;
        (-F*c_theta - m*L*omega^2*s_theta*c_theta - (M + m)*g*s_theta) / (L * denom)
    ];
end


% function sym_vector = createSymbolicVector(vecSymbol, n)
%     % Generate variable names dynamically
%     varNames = arrayfun(@(i) [vecSymbol, num2str(i)], 1:n, 'UniformOutput', false);
% 
%     % Define these as symbolic variables
%     syms(varNames{:});
% 
%     % Construct a symbolic vector from these variables
%     sym_vector = sym(varNames, 'real');
% 
%     sym_vector = sym_vector'; %column matrix by convention -- nx1
% end
% 
% %an infinite horizon LQR to stabilize the closed loop system at the terminal state
% function [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom_f, u_nom_f, Q, R)
% 
%     A_symb = jacobian(f, x);
%     B_symb = jacobian(f, u);
% 
%     %Get the linearised matrices (Jacobians) at terminal state
%     A_f = double(subs(A_symb, [x; u], [x_nom_f; u_nom_f]));
%     B_f = double(subs(B_symb, [x; u], [x_nom_f; u_nom_f]));
% 
%     [K_f,S_f,CLP] = lqr(A_f,B_f,Q,R)
% end


function plotPhasePortraits(trajectories, intialStatePerturbation, stateDims)
    % Plot all trajectories and nominal trajectory
    figure; hold on; grid on; axis equal;
    
    x0 = zeros(4,1);

    if nargin < 2
        stateDims = [1 3];
    end
    
    for i = 1:length(trajectories)
        x_sol = trajectories{i};
        plot(x_sol(stateDims(1), :), x_sol(stateDims(2), :), 'b-', 'LineWidth', 1.5);
    end

    xlabel(['x_{', num2str(stateDims(1)), '}'])
    ylabel(['x_{', num2str(stateDims(2)), '}'])
    title('MC Rollout Trajectories');
    
    drawnow

    %Visualise the set from which initial states were sampled
    % Generate angles for the circle
    theta = linspace(0, 2*pi, 100); % 100 points for a smooth circle
    
    % Calculate x and y coordinates
    x = intialStatePerturbation * cos(theta) + x0(stateDims(1));
    y = intialStatePerturbation * sin(theta) + x0(stateDims(2));
    
    % Plot the circle
    plot(x, y, '--g', 'LineWidth',1.4);
    

    %legend('Nominal Trajectory','IVP solution', 'Location','best');
    %hold off;
end


% Animated visualization function
function animate_cartpole_trajectory(t_opt, x_opt, params, saveVideoFlag, filename)
    % Animate the cartpole swing-up trajectory
    
    if nargin < 4
        saveVideoFlag = 0;
    end
    
    % if nargin < 5
    %     filename = 'cartpole_animation.mp4'; % default filename
    % end
    
    L = params.L;  % Pendulum length
    
    % Extract trajectories
    cart_x = x_opt(1,:);
    theta = pi - x_opt(3,:);
    
    % Calculate pole tip positions
    pole_tip_x = cart_x + L * sin(theta);
    pole_tip_y = L * cos(theta);
    
    % Set up figure
    figure('Position', [200, 200, 800, 600]);
    
    % Animation parameters
    playback_speed = 1.0;  % 1.0 = real time, 0.5 = half speed, 2.0 = double speed
    dt_anim = 0.05;        % Animation time step (seconds)
    
    % Time vector for animation
    t_anim = 0:dt_anim:(t_opt(end)/playback_speed);
    
    % Interpolate trajectories for smooth animation
    cart_x_anim = interp1(t_opt, cart_x, t_anim * playback_speed);
    theta_anim = interp1(t_opt, theta, t_anim * playback_speed);
    pole_tip_x_anim = cart_x_anim + L * sin(theta_anim);
    pole_tip_y_anim = L * cos(theta_anim);
    
    % Plot setup
    x_range = [min([cart_x, pole_tip_x]) - 0.3, max([cart_x, pole_tip_x]) + 0.3];
    y_range = [-L*1.2, L*1.2];
    
    % Initialize plot elements
    hold on; grid on; axis equal;
    xlim(x_range); ylim(y_range);
    xlabel('X Position (m)'); ylabel('Y Position (m)');
    title('Cartpole Swing-Up Animation');
    
    % Plot ground line
    plot(x_range, [0, 0], 'k-', 'LineWidth', 2);
    
    % Initialize animated elements
    h_cart = rectangle('Position', [0, -0.05, 0.2, 0.1], 'FaceColor', 'blue', 'EdgeColor', 'k');
    h_pole = plot([0, 0], [0, -L], 'r-', 'LineWidth', 4);
    h_tip = plot(0, -L, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    h_trail = plot(0, -L, 'b--', 'LineWidth', 1); %'MaxHeadSize', 0.5
    
    % Trail arrays
    trail_x = [];
    trail_y = [];
    max_trail_length = 50;  % Maximum number of trail points
    
    % % ==== Video writer setup ====
    % if saveVideoFlag
    %     v = VideoWriter(filename, 'MPEG-4');
    %     v.FrameRate = 1/dt_anim;   % match animation speed
    %     open(v);
    % end

    % Animation loop
    for i = 1:length(t_anim)
        % Current positions
        cart_pos = cart_x_anim(i);
        tip_x = pole_tip_x_anim(i);
        tip_y = pole_tip_y_anim(i);
        
        % Update cart position
        set(h_cart, 'Position', [cart_pos - 0.1, -0.05, 0.2, 0.1]);
        
        % Update pole
        set(h_pole, 'XData', [cart_pos, tip_x], 'YData', [0, tip_y]);
        
        % Update pole tip
        set(h_tip, 'XData', tip_x, 'YData', tip_y);
        
        % Update trail
        trail_x = [trail_x, tip_x];
        trail_y = [trail_y, tip_y];
        
        % Limit trail length
        if length(trail_x) > max_trail_length
            trail_x = trail_x(end-max_trail_length+1:end);
            trail_y = trail_y(end-max_trail_length+1:end);
        end
        
        set(h_trail, 'XData', trail_x, 'YData', trail_y);
        
        % Update title with current time and angle
        title(sprintf('Cartpole Swing-Up Animation (t=%.2fs, θ=%.1f°)', ...
              t_anim(i) * playback_speed, rad2deg(theta_anim(i))));
        
        % if saveVideoFlag
        %     % Capture frame and write to video
        %     frame = getframe(fig);
        %     writeVideo(v, frame);
        % end

        % Pause for animation
        pause(dt_anim);
        
        % Break if figure is closed
        if ~ishandle(h_cart)
            break;
        end
    end
    
    % Final message
    %if ishandle(h_cart)
    %    title('Cartpole Swing-Up Complete!');
    %end
    
    % if saveVideoFlag
    %     % Close video
    %     close(v);
    %     disp(['Animation saved as ', filename]);
    % end

end

function drawCircle(center, radius)
    %Visualise the set from which initial states were sampled
    % Generate angles for the circle
    theta = linspace(0, 2*pi, 100); % 100 points for a smooth circle
    
    % Calculate x and y coordinates
    x = radius * cos(theta) + center(1);
    y = radius * sin(theta) + center(2);
    
    % Plot the circle
    plot(x, y, '--g', 'LineWidth',1.4);
end