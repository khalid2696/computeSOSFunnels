%clc; clearvars; close all
<<<<<<< HEAD
%yalmip('clear');
=======
yalmip('clear');
>>>>>>> unicycle

disp('Hang on..');
disp('Computing nominal trajectory using direct collocation trajectory optimisation');
disp(' ');

%% Add directories
addpath('./lib/');

%% Inherit parameters if they exist in the wrapper file

% Initial and target states : [pos, vel, theta, omega]
if exist('initialState','var') || exist('finalState','var')
    x0 = initialState;
    xf = finalState;
else
<<<<<<< HEAD
    x0 = [0; 0; 0; 0;];   %Initial state: pendulum hanging down, cart at origin
    xf = [1; 0; pi; 0;];  %Final state: pendulum balanced upright, cart at specific position
=======
    x0 = [0; 0; pi/2;];    %random initial state
    xf = [0; 2; pi/2; 0;]; %random final state
>>>>>>> unicycle
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
    drawFlag = 1; %by default do not plot
end

%% Nominal trajectory and nominal/open-loop control input (feedforward term)

tic
[x_nom, u_nom, time_instances, nom_trajCost, diagnostics] = getNominalTrajectory_using_DirectCollocation(dynamicsFnHandle, x0, xf, T_max, N);
toc   

disp('Finshed computing nominal trajectory and nominal (feedforward) input tape');
disp(' ');

if diagnostics.problem ~= 0 %if the optimisation fails, relax the time-horizon and input constraints
    disp(diagnostics.info);
    error('Failed to compute a feasible nominal trajectory! Cannot proceed.. Exitting!');
end

%% save the nominal trajectory and feedforward input
save('./precomputedData/nominalTrajectory.mat', 'time_instances', 'x_nom', 'u_nom', 'nom_trajCost', 'dynamicsFnHandle');

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

    plotStateInputTrajectories(time_instances, x_nom, u_nom)
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

function plotStateInputTrajectories(t_opt, x_opt, u_opt)
    figure;
        
    subplot(2,2,1);
    plot(x_opt(1,:), x_opt(2,:), 'b-', 'LineWidth', 2);
    xlabel('p_x [m]'); ylabel('p_y [m]');
    title('Unicycle position'); grid on;
    
    subplot(2,2,2);
    plot(t_opt, rad2deg(x_opt(3,:)), 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel(' \theta [deg]');
    title('Unicycle orientation vs Time'); grid on;
    
    subplot(2,2,3);
    plot(t_opt, u_opt(1,:), 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel(' v [m/s]');
    title('Input Velocity'); grid on;
    
    subplot(2,2,4);
    plot(t_opt, u_opt(2,:), 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel(' \omega [rad/s]');
    title('Input angular velocity'); grid on;

end
