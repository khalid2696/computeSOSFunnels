%clc; clearvars; close all
%yalmip('clear');

disp('Hang on..');
disp('Computing nominal trajectory using direct collocation trajectory optimisation');
disp(' ');

% Quadrotor behavior:
% -ve roll for +ve y
% +ve pitch for +ve x
% zB vertically up

%% Add directories
addpath('./lib/');

%% Inherit parameters if they exist in the wrapper file

% if ~exist('quadParameters','var')
%     quadParameters.m = 0.5;        % mass (kg)
%     quadParameters.g = 9.81;       % gravity (m/s^2)
%     quadParameters.J = [5e-3, 5e-3, 10e-3]; % moment of inertia (kg⋅m^2)                    
% end

% Initial & target pose: position and Euler angles (roll, pitch, yaw) in radians
if exist('initialPose','var') || exist('finalPose','var')
    p0 = initialPose(1:3); euler0 = initialPose(4:6);
    pf = finalPose(1:3); eulerf = finalPose(4:6);
else
    p0 = [0; 0; 2]; euler0 = [0; 0; 0];  % default initial pose -- at origin pointing North 
    pf = [1.5; 4.5; 2]; eulerf = [0; 0; 0];  % some random final pose
end

if exist('numTimeSteps','var')
    N = numTimeSteps;
else
    N = 25; % default number of time steps
end

% Maximum allowable time horizon
if exist('maxTimeHorizon','var')
    T_max = maxTimeHorizon;          
else
    T_max = 10;  %some default value
end

if ~exist('drawFlag','var')
    drawFlag = 0; %by default do not plot
end

%% Initial and target states : [pos, vel, angles, omega]

x0 = [p0; zeros(3,1); euler0; zeros(3,1)]; %start from zero linear/angular velocities
xf = [pf; zeros(3,1); eulerf; zeros(3,1)]; %end at zero linear/angular velocities

%% Nominal trajectory and nominal/open-loop control input (feedforward term)

tic
[x_nom, u_nom, time_instances, nom_trajCost, diagnostics] = ...
        getNominalTrajectory_using_DirectCollocation(dynamicsFnHandle, x0, xf, T_max, N, quadParameters);
toc

disp('Finshed computing nominal trajectory and nominal (feedforward) input tape');
disp(' ');

if diagnostics.problem ~= 0 %if the optimisation fails, relax the time-horizon and input constraints
    disp(diagnostics.info);
    error('Failed to compute a feasible nominal trajectory! Cannot proceed.. Exitting!');
end



%% save the nominal trajectory and feedforward input
save('./precomputedData/nominalTrajectory.mat', 'time_instances', 'x_nom', 'u_nom', 'dynamicsFnHandle', 'quadParameters');

disp('Saved the nominal trajectory and nominal inputs to a file!');
disp(' ');

%% Additionally, plot if enabled

if drawFlag
    
    figure; plot_flat_outputs(time_instances, x_nom, 2.5); %last argument: scaling for the quadrotor visual
    
    % ---- Just for dev and debug -----%
    figure; plot_Euler_angles(time_instances, x_nom);
    figure; plot_linear_velocities(time_instances, x_nom)
    figure; plot_angular_velocities(time_instances, x_nom)

    figure; plot_input_profiles(time_instances, u_nom);
end

clearvars; %cleanup the workspace after saving relevant data and plotting

% %% for debugging purposes
% drawFlag = 0;
% 
% load('./precomputedData/nominalTrajectory.mat');
% 
% tspan = [time_instances(1), time_instances(end)];
% 
% % Solve IVP
% x0 = x_nom(:,1);
% [t_sol, x_sol] = ode15s(@(t,x) ...
%     ode_wrapper(t, x, time_instances, u_nom', dynamicsFnHandle), ...
%     tspan, x0);
% t_sol = t_sol';
% x_sol = x_sol';
% u_sol = interp1(time_instances, u_nom', t_sol, 'linear')';
% 
% x_coarse = interp1(t_sol, x_sol', time_instances, 'spline')';
% 
% if drawFlag
% 
%     figure; plot_flat_outputs(time_instances, x_nom, 2.5); %last argument: scaling for the quadrotor visual
%     figure; plot_flat_outputs(t_sol, x_sol, 2.5);
% 
%     % ---- Just for dev and debug -----%
%     figure; plot_Euler_angles(time_instances, x_nom);
%     figure; plot_Euler_angles(t_sol, x_sol);
% 
%     figure; plot_linear_velocities(time_instances, x_nom);
%     figure; plot_linear_velocities(t_sol, x_sol);
% 
%     figure; plot_angular_velocities(time_instances, x_nom)
%     figure; plot_angular_velocities(t_sol, x_sol);
% 
%     figure; plot_input_profiles(time_instances, u_nom);
%     figure; plot_input_profiles(t_sol, u_sol);
% end
% 
% [x_fine, u_fine] = upsample_state_control_trajectories(time_instances, x_nom, u_nom, t_sol);

%% Local Function definitions

% function q = eul2quat(eulerAngle)
%     % Convert Euler angles [roll; pitch; yaw] to quaternion [w; x; y; z]
%     roll = eulerAngle(1); pitch = eulerAngle(2); yaw = eulerAngle(3);
% 
%     cr = cos(roll/2); sr = sin(roll/2);
%     cp = cos(pitch/2); sp = sin(pitch/2);
%     cy = cos(yaw/2); sy = sin(yaw/2);
% s
%     q = [cr*cp*cy + sr*sp*sy; % w
%          sr*cp*cy - cr*sp*sy; % x
%          cr*sp*cy + sr*cp*sy; % y
%          cr*cp*sy - sr*sp*cy]; % z
% end
% 
% function euler = quat2euler(q)
%     % Converts quaternion [q0, q1, q2, q3] to Euler angles [roll; pitch; yaw]
%     q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
% 
%     % Compute Euler angles
%     roll = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1^2 + q2^2)); % φ (roll)
%     pitch = asin(2 * (q0 * q2 - q3 * q1)); % θ (pitch)
%     yaw = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2^2 + q3^2)); % ψ (yaw)
% 
%     % Return as a column vector
%     euler = [roll; pitch; yaw];
% end

% ODE wrapper function
function xdot = ode_wrapper(t, x, time_instances, u_nom, dynamicsFnHandle)
    
    u = interp1(time_instances, u_nom, t, 'previous', 'extrap');  % hold last value
    
    xdot = dynamicsFnHandle(x, u);
end

%Interpolates state and control vectors
function [x_fine, u_fine, t_fine] = upsample_state_control_trajectories(t, x, u, t_fine)

    % Dimensions
    Nd = length(t_fine); n = size(x, 1); m = size(u, 1);

    % Preallocate outputs
    x_fine = zeros(n, Nd);
    u_fine = zeros(m, Nd);

    % Interpolate each row (dimension) of x with cubic spline
    for i = 1:n
        x_fine(i, :) = interp1(t, x(i, :), t_fine, 'spline'); % Cubic interpolation
    end

    % Interpolate each row (dimension) of u with linear interpolation
    for i = 1:m
        u_fine(i, :) = interp1(t, u(i, :), t_fine, 'linear'); % Linear interpolation
    end
end


function plot_flat_outputs(time_instances, x_nom, scaling)

    if nargin < 3
        scaling = 1; %scaling for the quadrotor visual
    end
    
    %Plot results
    subplot(3,1,1)
    plot(x_nom(1,:),x_nom(2,:),'--k'); %x-y location
    axis equal; hold on;

    xlabel('p_x [m]');
    ylabel('p_y [m]');

    %initial pose
    plot_quadrotor(x_nom(1,1), x_nom(2,1), pi/2, scaling); %x,y, heading   
    % Add a marker for the center position
    plot(x_nom(1,1), x_nom(2,1), 'sg', 'MarkerSize', 3, 'MarkerFaceColor', 'g');

    %final pose
    plot_quadrotor(x_nom(1,end), x_nom(2,end), pi/2, scaling); %x,y, heading   
    % Add a marker for the center position
    plot(x_nom(1,end), x_nom(2,end), 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
    title('Quadrotor Position and Heading');

    subplot(3,1,2)
    plot(time_instances, x_nom(3,:)); %altitude
    xlabel('time [s]');
    ylabel('p_z [m]');
    
    subplot(3,1,3)
    plot(time_instances, x_nom(9,:)*180/pi); %heading
    xlabel('time');
    ylabel('Heading [deg]');
end


function plot_quadrotor(x, y, heading, scaling)
    
    if nargin < 4
        scaling = 1;
    end

    % Parameters
    arm_length = 0.1*scaling;
    arrow_length = 0.25*scaling;
    rotor_radius = 0.04*scaling;

    % Define arm endpoints relative to center (cross)
    arms = [ -arm_length,  0;
              arm_length,  0;
              0, -arm_length;
              0,  arm_length ]';

    % Rotation matrix for heading
    R = [cos(heading), -sin(heading); sin(heading), cos(heading)];

    % Rotate arms
    rotated_arms = R * arms;

    % Translate to (x, y)
    quad_x = rotated_arms(1,:) + x;
    quad_y = rotated_arms(2,:) + y;

    % Plot quadrotor arms (horizontal and vertical)
    plot(quad_x(1:2), quad_y(1:2), 'b-', 'LineWidth', 1.5); hold on;
    plot(quad_x(3:4), quad_y(3:4), 'b-', 'LineWidth', 1.5);

    % Plot heading arrow
    quiver(x, y, arrow_length*cos(heading), arrow_length*sin(heading), 0, ...
           'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);

    % Plot center point
    plot(x, y, 'ko', 'MarkerFaceColor', 'k');

    % Plot rotors as circles at the ends of the arms
    theta = linspace(0, 2*pi, 50);
    for i = 1:4
        rotor_x = rotor_radius * cos(theta) + quad_x(i);
        rotor_y = rotor_radius * sin(theta) + quad_y(i);
        fill(rotor_x, rotor_y, 'c', 'EdgeColor', 'k'); % Green rotors with black edge
    end
end


function plot_input_profiles(time_instances_opt, u_opt)

    subplot(4,1,1)
    plot(time_instances_opt, u_opt(1,:)); %total thrust
    ylabel('Thrust [N]');
    
    subplot(4,1,2)
    plot(time_instances_opt, u_opt(2,:)); %Roll moment
    ylabel('M_p [N-m]');
    
    subplot(4,1,3)
    plot(time_instances_opt, u_opt(3,:)); %Pitch moment
    ylabel('M_q [N-m]');
    
    subplot(4,1,4)
    plot(time_instances_opt, u_opt(4,:)); %Yaw moment
    ylabel('M_r [N-m]');
end

function plot_Euler_angles(time_instances_opt, x_opt)
    

    % Euler angles
    title('Euler angles');
    
    subplot(3,1,1)
    plot(time_instances_opt, x_opt(7,:)*180/pi); %phi
    xlabel('time');
    ylabel('phi [deg]');
    
    subplot(3,1,2)
    plot(time_instances_opt, x_opt(8,:)*180/pi); %theta
    xlabel('time');
    ylabel('theta [deg]');
    
    subplot(3,1,3)
    plot(time_instances_opt, x_opt(9,:)*180/pi); %psi
    xlabel('time');
    ylabel('psi [deg]');
end

function plot_linear_velocities(time_instances_opt, x_opt)
    

    % Euler angles
    title('Linear Velocities');
    
    subplot(3,1,1)
    plot(time_instances_opt, x_opt(4,:)); %v_x
    xlabel('time');
    ylabel('v_x [m/s]');
    
    subplot(3,1,2)
    plot(time_instances_opt, x_opt(5,:)); %v_y
    xlabel('time');
    ylabel('v_y [m/s]');
    
    subplot(3,1,3)
    plot(time_instances_opt, x_opt(6,:)); %v_z
    xlabel('time');
    ylabel('v_z [m/s]');
end

function plot_angular_velocities(time_instances_opt, x_opt)
    

    % Euler angles
    title('Angular velocities');
    
    subplot(3,1,1)
    plot(time_instances_opt, x_opt(10,:)); %p
    xlabel('time');
    ylabel('p [rad/s]');
    
    subplot(3,1,2)
    plot(time_instances_opt, x_opt(11,:)); %q
    xlabel('time');
    ylabel('q [rad/s]');
    
    subplot(3,1,3)
    plot(time_instances_opt, x_opt(12,:)); %r
    xlabel('time');
    ylabel('r [rad/s]');
end

