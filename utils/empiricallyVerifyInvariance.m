%clc; clearvars; 
%close all

%% Add directories
addpath('../lib/');

%% Load the precomputed files

load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/LQRGainsAndCostMatrices.mat');
load('../precomputedData/setInvarianceCertificates.mat')

%% Status message
%Quantities at our disposal now

% N                  : number of time samples                   : scalar N
% time_instances     : time horizon (sampled)                   : 1 x N
% x_nom              : nominal state trajectory                 : n_x x N
% u_nom              : feedforward input tape                   : n_u x N
% K                  : Feedback gains (sampled)                 : n_u x n_x x N
% P                  : cost-to-goal matrix (sampled)            : n_x x n_x x N
% candidateV         : Lyapunov certificates of invariance      : 1 x N (cell)
% ellipsoidMatrices  : ellipsoids characterizing candidateV     : n_x x n_x x N 
% rhoScaling         : level-set boundary value                 : 1 x N
% terminalRegion     : Ellipsoidal goal region in BRS analysis  : n x n
% multiplierTerms    : polynomial multipliers from S-procedure  : 1 x N-1 (cell)

N = length(time_instances);
n = size(x_nom, 1); m = size(u_nom, 1); %state and input vector dimensionality

%% Specify parameters or Inherit them if they exist in the wrapper file

%Specify the start time for the rollouts
if ~exist('startTimeIndex','var')
    startTimeIndex = 1;
end

if ~exist('numSamples','var')
    numSamples = 500; %default number of rollouts
end

if ~exist('projectionDims_2D','var')
    projectionDims_2D = [1 3];
end 

%% Monte Carlo rollouts for empirical verification

%t0 = time_instances(1); tf = time_instances(end); 
t_start = time_instances(startTimeIndex);
Ts = (time_instances(end) - time_instances(1))/(N-1); %sampling time

rollout_time_horizon = time_instances(startTimeIndex:end);
rollout_x_nom        = x_nom(:,startTimeIndex:end);
rollout_u_nom        = u_nom(:,startTimeIndex:end);
rollout_K            = K(:,:,startTimeIndex:end);
rollout_P            = P(:,:,startTimeIndex:end);
rollout_ROA          = ellipsoidMatrices(:,:,startTimeIndex:end);
rollout_levelSetVals = rhoScaling(startTimeIndex:end);

%sample initial states at random
intialState = rollout_x_nom(:,1); %centered around the nominal trajectory
inletSet = ellipsoidMatrices(:,:,startTimeIndex)/rhoScaling(startTimeIndex);

%add argument: 'surface' to the following function for sampling only on the ellipsoid surface
initialStateSamples = samplePointsFromEllipsoid(intialState, inletSet, numSamples);

trajectories = cell(numSamples, 1);
input_profiles = cell(numSamples, 1);
LyapunovFnValue_profiles = cell(numSamples, 1);

errors = zeros(numSamples, length(rollout_time_horizon));
costs = zeros(numSamples, length(rollout_time_horizon));

for i = 1:numSamples
    x0 = initialStateSamples(:, i);
    %options: Euler, trapezoidal, RK4, (inbuilt) ode45
    [x_traj, total_input, error, cost] = ...
        forward_propagate(dynamicsFnHandle, x0, rollout_x_nom, rollout_u_nom, rollout_K, rollout_P, rollout_time_horizon, Ts, 'RK4');
    trajectories{i} = x_traj;
    input_profiles{i} = total_input;
    errors(i, :) = error;
    costs(i, :) = cost;

    %computing Lyapunov function value
    LyapunovFnValsArray = NaN(size(rollout_time_horizon));
    for k=1:length(rollout_time_horizon)
        tempCurrState = rollout_x_nom(:,k);
        tempEllipsoidMatrix = rollout_ROA(:,:,k)/rollout_levelSetVals(k);
        LyapunovFnValsArray(k) = (x_traj(:,k)-tempCurrState)'*tempEllipsoidMatrix*(x_traj(:,k)-tempCurrState);
    end

    LyapunovFnValue_profiles{i} = LyapunovFnValsArray;
end

disp('-- End of Monte Carlo empirical verification --');
disp(' ');

%% Visualization
disp('Plotting trajectories and metrics from MC rollouts..');
disp(' ');

%Plot the rollout trajectories
plot_trajectories(trajectories, rollout_x_nom, inletSet, x_nom, projectionDims_2D);

%Plot the input profiles
%plot_input_profiles(input_profiles, rollout_time_horizon, u_nom, time_instances)

%Visualise metrics
plotMetrics(LyapunovFnValue_profiles, errors, costs, rollout_time_horizon);

clearvars;

%% Function defintions

% Sampling interior or surface of an Ellipsoid   
function samplePoints = samplePointsFromEllipsoid(center, ellipsoidMatrix, numSamplePoints, samplingRegion)
    
    if nargin < 4
        samplingRegion = 'interior';
    end
    
    % Eigen decomposition
    [Q, Lambda] = eig(ellipsoidMatrix);
    n = size(center,1);

    % Semi-axis lengths
    semi_axes_lengths = 1 ./ sqrt(diag(Lambda));
    
    % Generate random points on the n-dimensional unit sphere
    sphere_points = randn(n, numSamplePoints); % Random points
    sphere_points = sphere_points ./ vecnorm(sphere_points); % Normalize to lie on the unit sphere
    
    if strcmpi(samplingRegion,'interior')
        scaling = rand(n,numSamplePoints); %multiply by scaling between 0 and 1 for interior
    else
        scaling = ones(n,numSamplePoints); %multiply by scaling of 1 for surface
    end

    sphere_points = sphere_points .* scaling;

    % Transform points to the ellipsoid
    samplePoints = Q * diag(semi_axes_lengths) * sphere_points + center;
end

function [x_traj, total_input, errorNorm, costToGoal] = forward_propagate(dynamics, x0, x_nom, u_nom, K, P, time, dt, method)
    % Simulates the unicycle dynamics under TVLQR control
    x_traj = zeros(size(x_nom));
    x_traj(:, 1) = x0;
    total_input = zeros(size(u_nom));

    errorNorm = zeros(1, length(time));
    costToGoal = zeros(1, length(time));
    
    for k = 1:length(time)-1
        u_k = u_nom(:, k) - K(:, :, k) * (x_traj(:, k) - x_nom(:, k));
        
        % Propagate dynamics
        if strcmpi(method,'Euler')
            %Simpler Euler integration -- for speed
            dx = dynamics(x_traj(:, k), u_k);
            x_traj(:, k+1) = x_traj(:, k) + dt * dx;
        elseif strcmpi(method,'trapezoidal')
            % Trapezoidal integration -- for accuracy and speed
            dx1 = dynamics(x_traj(:, k), u_k);               % Slope at start of interval
            x_temp = x_traj(:, k) + dt * dx1;                         % Euler step for estimation
            dx2 = dynamics(x_temp, u_k);                     % Slope at end of interval
            x_traj(:, k+1) = x_traj(:, k) + (dt / 2) * (dx1 + dx2);   % Trapezoidal integration
        elseif strcmpi(method,'RK4')
            %RK4 integration
            k1 = dynamics(x_traj(:, k), u_k); 
            k2 = dynamics(x_traj(:, k) + 0.5 * dt * k1, u_k);
            k3 = dynamics(x_traj(:, k) + 0.5 * dt * k2, u_k);
            k4 = dynamics(x_traj(:, k) + dt * k3, u_k);
            x_traj(:, k+1) = x_traj(:, k) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4); % Update state
        else %just use the in-built ode45 (slower but accurate)
            [~, x_next] = ode45(@(t, x) dynamics(x, u_k), [0, dt], x_traj(:, k));
            x_traj(:, k+1) = x_next(end, :)';
        end
    
        % Compute error and cost metrics
        stateDeviation = x_traj(:, k) - x_nom(:, k);
        errorNorm(k) = norm(stateDeviation);
        costToGoal(k) = stateDeviation' * P(:, :, k) * stateDeviation;
        total_input(:,k) = u_k;
    end
    
    total_input(:,end) = u_nom(:,end);
    finalStateDeviation = x_traj(:, end) - x_nom(:, end);
    errorNorm(end) = norm(finalStateDeviation);
    costToGoal(end) = finalStateDeviation' * P(:, :, end) * finalStateDeviation;
end


function plot_trajectories(trajectories, rollout_x_nom, complete_x_nom, projectionDims)
    
    if nargin < 4
        projectionDims = [1 2]; %if not specified, by default x-y projection
    end

    % Plot all trajectories and nominal trajectory
    %figure; hold on; grid on; axis equal;
    
    plot(rollout_x_nom(projectionDims(1), :), rollout_x_nom(projectionDims(2), :), 'k--', 'LineWidth', 2);

    for i = 1:length(trajectories)
        plot(trajectories{i}(projectionDims(1), :), trajectories{i}(projectionDims(2), :), 'b-', 'LineWidth', 0.5);
        plot(trajectories{i}(projectionDims(1), 1), trajectories{i}(projectionDims(2), 1), 'sg');
    end
    
    plot(complete_x_nom(projectionDims(1), :), complete_x_nom(projectionDims(2), :), 'k--', 'LineWidth', 2);

    xlabel(['x_{', num2str(projectionDims(1)), '}']);
    ylabel(['x_{', num2str(projectionDims(2)), '}']);
    title('Monte Carlo Rollout Trajectories');
    %legend('Nominal Trajectory','Trajectory rollouts','Sample Initial states','Location','best');
end

function plot_input_profiles(input_profiles, rollout_time_instances, complete_u_nom, complete_time_instances)
    % Plot all trajectories and nominal trajectory
    figure;
    
    subplot(2,1,1)
    hold on;
    grid on; 
    for i = 1:length(input_profiles)
        plot(rollout_time_instances, input_profiles{i}(1, :), 'b-', 'LineWidth', 0.5);
    end
    plot(complete_time_instances, complete_u_nom(1, :), 'k--', 'LineWidth', 2);
    xlabel('time'); ylabel('u_1');
    title('Input history from Monte Carlo Rollouts');
    
    subplot(2,1,2)
    hold on;
    grid on; 
    for i = 1:length(input_profiles)
        plot(rollout_time_instances, input_profiles{i}(2, :), 'b-', 'LineWidth', 0.5);
    end
    plot(complete_time_instances, complete_u_nom(2, :), 'k--', 'LineWidth', 2);

    xlabel('time'); ylabel('u_2');
    %legend('Sample Trajectories', 'Initial states', 'Nominal Trajectory', '');
end

% Plot error and cost metrics over time
function plotMetrics(LyapunovFnValue_profiles, errors, costs, time)
    
    figure;
     
    %Lyapunov function value -- 1-sublevel-set value
    subplot(3,1,1);
    hold on

    for i=1:length(LyapunovFnValue_profiles)
        thisProfile = LyapunovFnValue_profiles{i};
        plot(time,thisProfile,'b-','LineWidth', 1.2);
    end
    
    plot(time, ones(size(time)), '--k', 'LineWidth', 1.5)
    grid on
    xlim([0 time(end)+0.5]); ylim([0 1.2]);
    xlabel('Time (s)'); ylabel('Function Value');
    title('Lyapunov function value of 1-sublevel-set');
    

    %error-norm
    subplot(3, 1, 2);
    hold on
    
    % Calculate upper and lower bounds of data spread
    upper_bound = mean(errors, 1) + std(errors,1);
    lower_bound = mean(errors, 1) - std(errors,1);

    % Plot the shaded std dev region
    fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
        [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Shaded region
    
    %Plot the mean
    plot(time, mean(errors, 1), 'r-', 'LineWidth', 2);
    grid on;
    xlim([0 time(end)+0.5]);
    xlabel('Time (s)'); ylabel('Error Norm');
    title('Mean of Error Norm Over Time');
    

    %cost-to-goal metric (weighted with P_k)
    subplot(3, 1, 3);
    hold on
    
    % Calculate upper and lower bounds of data spread
    upper_bound = mean(costs, 1) + std(costs,1);
    lower_bound = mean(costs, 1) - std(costs,1);
    % Plot the shaded std dev region
    fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
        [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Shaded region
    
    %Plot the mean
    plot(time, mean(costs, 1), 'b-', 'LineWidth', 2);
    grid on;
    xlim([0 time(end)+0.5]);
    xlabel('Time (s)'); ylabel('Weighted Cost');
    title('Mean of Cost-to-Goal Metric Over Time');
end


function plotTerminalSet(x_final, terminalEllipsoid)

    M = terminalEllipsoid;
    M_xy = project_ellipsoid_matrix(M, [1 2]);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_final(1) + rotated_ellipse(1, :), ...
         x_final(2) + rotated_ellipse(2, :), ...
         '-.r', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
end

function M_2D = project_ellipsoid_matrix(M, projection_dims)
    % Input:
    % M: nxn matrix defining the n-dimensional ellipsoid x^T M x < 1
    % projection_dims: 2-element vector specifying which dimensions to project onto
    %                  (e.g., [1 2] for xy-plane, [1 3] for xz-plane)
    
    n = size(M, 1); %get the dimensionality of matrix M

    basisMatrix = zeros(n,2);
    basisMatrix(projection_dims(1),1) = 1;
    basisMatrix(projection_dims(2),2) = 1;

    M_2D = inv(basisMatrix' *inv(M) * basisMatrix);
end

function plotEllipse(center, ellipseMatrix)
    
    %plot an ellipse from which initial states are sampled
    ellipseCenter = center(1:2); % 2D center of the ellipsoid
    [eig_vec, eig_val] = eig(ellipseMatrix);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(ellipseCenter(1) + rotated_ellipse(1, :), ...
         ellipseCenter(2) + rotated_ellipse(2, :), ...
         '-k', 'LineWidth', 1.2);  
end