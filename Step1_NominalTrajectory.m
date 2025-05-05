%clc; clearvars; close all

%% Add directories
addpath('./lib/');

%% Inherit parameters if they exist in the wrapper file

if exist('numTimeSteps','var')
    N = numTimeSteps;
else
    N = 25; % default number of time steps
end

if ~exist('drawFlag','var')
    drawFlag = 0; %by default do not plot
end

% Maximum allowable time horizon
if exist('maxTimeHorizon','var')
    T_max = maxTimeHorizon;          
else
    T_max = 10;  %some default value
end

%Initial and Final states for constraints: [x; y; theta]
if exist('xinitial','var') || exist('xfinal','var')
    x0 = xinitial;
    xf = xfinal;
else
    x0 = [0; 0; pi/2];   % default initial state -- at origin pointing North 
    xf = [2; 5; pi/2];   % some random final state
end
%% System dynamics - define states, input and nonlinear function
% Unicycle dynamics: x_dot = f(x, u)
dynamicsFnHandle = @(x, u) [u(1) * cos(x(3)); ...  % x_dot
                            u(1) * sin(x(3)); ...  % y_dot
                             u(2)];                % theta_dot

%% Nominal trajectory and nominal/open-loop control input (feedforward term)

% Trajectory Optimisation Parameters
%N = 25;              % Number of time steps
%T_max = 10;          % Maximum allowable time horizon
input_saturations = [0.75; pi/2]; %Linear and angular velocity limits (absolute value)
                             % -input_saturations(j) <= u_j <= input_saturations(j) 
%input_saturations = [1; pi]; %Linear and angular velocity limits (absolute value)
                             % -input_saturations(j) <= u_j <= input_saturations(j) 

[x_nom, u_nom, time_instances, nom_trajCost, diagnostics] = ...
getNominalTrajectory_using_DirectCollocation(dynamicsFnHandle, x0, xf, input_saturations, T_max, N);

disp('Finshed computing nominal trajectory and nominal (feedforward) input tape');
disp(' ');

if diagnostics.problem ~= 0 %if the optimisation fails, relax the time-horizon and input constraints
    disp(diagnostics.info);
    error('Failed to compute a feasible nominal trajectory! Cannot proceed.. Exitting!');
end

%% save the nominal trajectory and feedforward input
save('./precomputedData/nominalTrajectory.mat', 'time_instances', 'x_nom', 'u_nom', 'dynamicsFnHandle');

disp('Saved the nominal trajectory and nominal inputs to a file!');
disp(' ');

%% Additionally, plot if enabled

if drawFlag
    figure; hold on; grid on; axis equal

    %x-y trajectory
    plot(x_nom(1,:),x_nom(2,:),'-.k','LineWidth',1.2);     
    
    %initial pose
    plot_pentagon_car(x_nom(:,1));   
    % Add a marker for the center position
    plot(x_nom(1,1), x_nom(2,1), 'sg', 'MarkerSize', 3, 'MarkerFaceColor', 'g');

    %final pose
    plot_pentagon_car(x_nom(:,end));   
    % Add a marker for the center position
    plot(x_nom(1,end), x_nom(2,end), 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r');

    xlabel('p_x'); ylabel('p_y');
end

clearvars; %cleanup the workspace after saving relevant data and plotting

%% Local function definitions

% Function to plot a pentagon representing a car at position (x, y) with orientation theta.
function plot_pentagon_car(pose)
    
    x = pose(1); y = pose(2); theta = pose(3) - pi/2;
    % Pentagon vertices in local coordinates (centered at origin)
    pentagon_local = [
        0,  0.5;  % Front vertex
        -0.3, 0.2; % Front-left
        -0.3, -0.3; % Rear-left
        %0, -0.5;  % Rear
        0.3, -0.3; % Rear-right
        0.3, 0.2;  % Front-right (connect back to front-left)
    ]';

    % Rotation matrix for orientation theta
    R = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)];

    % Scale, Rotate and translate the pentagon
    scaling = 0.5;
    pentagon_world = scaling * R * pentagon_local + [x; y];

    % Plot the pentagon
    fill(pentagon_world(1, :), pentagon_world(2, :), 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    axis equal;
end