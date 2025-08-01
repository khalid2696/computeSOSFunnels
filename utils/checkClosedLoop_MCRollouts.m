%clc; clearvars; close all

%% Add directories
addpath('../lib/');

%% Load the nominal trajectory and LQR gains
load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/LQRGainsAndCostMatrices.mat');
N = length(time_instances);

% -- Status message: Quantities at our disposal now -- %

% time_instances   : time horizon (sampled)        : 1 x N
% x_nom            : nominal state trajectory      : n_x x N
% u_nom            : nominal input tape            : n_u x N
% K                : Feedback gains (sampled)      : n_u x n_x x N
% P                : cost-to-goal matrix (sampled) : n_x x n_x x N
% dynamicsFnHandle : the function handle of the system dynamics

%% Specify parameters or Inherit them if they exist in the wrapper file

%Specify the start time for the rollouts
if ~exist('startTimeIndex','var')
    %startTimeIndex = randi(N-1); %should be in between 1 and N-1
    startTimeIndex = 1;
end

if ~exist('num_samples','var')
    num_samples = 500; %default number of rollouts
end

if exist('startMaxPerturbation','var')
    rho0 = startMaxPerturbation;
else
    rho0 = 0.3; %decrease this for a smaller initial set
end

%finer discretization to prevent integration error build-up
if ~exist('upsamplingFactor','var')
    upsamplingFactor = 40;
end

%% Upsample trajectories and matrices for forward rollouts

Nd = N*upsamplingFactor; % number of interpolated points

% Fine time vector for interpolation: upsampled spacing
t_fine = linspace(time_instances(1), time_instances(end), Nd);

[x_nom, u_nom] = upsample_state_control_trajectories(time_instances, x_nom, u_nom, t_fine);

P = upsample_matrix(P, time_instances, t_fine);
K = upsample_matrix(K, time_instances, t_fine);

time_instances = t_fine;

%% Monte Carlo forward rollouts

% a scaling to vary initial set size based on the time it starts
% required because we do systems analysis over a *finite-time* horizon
t0 = time_instances(1); tf = time_instances(end); t_start = time_instances(startTimeIndex);
Ts = (time_instances(end)-time_instances(1))/(length(time_instances)-1); %sampling-time

c = 3;      %increase this for steeper decrease in rho(t)
rho = rho0*exp(c*(t_start-t0)/(t0-tf));

%Sampling initial states from an initial ellipsoidal set
if ~exist('initial_state_covariance', 'var')
    initial_state_covariance = (1/3)*rho*P(:,:,startTimeIndex)^(-1/2);
    %mu + 3*sigma will cover 99.7% of the distribution
    %hence divided by 3
end

%1/sqrt(eigenvalue of P_k) gives the length of each semi-axis
% deviation^T P_k deviation = rho
% and eigen value of P_k^(-1/2) is 1/sqrt(eigenvalue of P_k)
% that is why, we take P^(-1/2)

rollout_time_horizon = time_instances(startTimeIndex:end);
rollout_x_nom = x_nom(:,startTimeIndex:end);
rollout_u_nom = u_nom(:,startTimeIndex:end);
rollout_K     = K(:,:,startTimeIndex:end);
rollout_P     = P(:,:,startTimeIndex:end);

%sample initial states at random
mean_intial_state = rollout_x_nom(:,1); %centered around the nominal trajectory

initial_states = sample_initial_states(mean_intial_state, initial_state_covariance, num_samples);

trajectories = cell(num_samples, 1);
input_profiles = cell(num_samples, 1);

errors = zeros(num_samples, length(rollout_time_horizon));
costs = zeros(num_samples, length(rollout_time_horizon));

for i = 1:num_samples
    x0 = initial_states(:, i);
    %options: Euler, trapezoidal, RK4, (inbuilt) ode45
    [x_traj, total_input, error, cost] = ...
        forward_propagate(dynamicsFnHandle, x0, rollout_x_nom, rollout_u_nom, rollout_K, rollout_P, rollout_time_horizon, Ts,'RK4');
    trajectories{i} = x_traj;
    input_profiles{i} = total_input;
    errors(i, :) = error;
    costs(i, :) = cost;
end

disp('-- End of Monte Carlo forward rollouts --');
disp(' ');

%% Visualization
disp('Plotting trajectories, input profiles, and metrics from MC rollouts..');
disp(' ');

plot_xy_trajectories(trajectories, rollout_x_nom, initial_state_covariance, x_nom);
%plot_xyz_trajectories(trajectories, rollout_x_nom, initial_state_covariance, x_nom);

%plot state trajectories
plot_state_trajectories(trajectories, rollout_time_horizon, rollout_x_nom, [1 2 3]);    %position
plot_state_trajectories(trajectories, rollout_time_horizon, rollout_x_nom, [4 5 6]);    %velocity
plot_state_trajectories(trajectories, rollout_time_horizon, rollout_x_nom, [7 8 9]);    %attitude
plot_state_trajectories(trajectories, rollout_time_horizon, rollout_x_nom, [10 11 12]); %body rates

plot_input_profiles(input_profiles, rollout_time_horizon, u_nom, time_instances);
plot_error_metrics(errors, costs, rollout_time_horizon);

clearvars;

%% Function defintions

function initial_states = sample_initial_states(mean_state, covariance, num_samples)
    % Samples initial states from an ellipsoid
    initial_states = mvnrnd(mean_state, covariance, num_samples)';
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

%Interpolates state and control vectors
function [x_fine, u_fine, t_fine] = upsample_state_control_trajectories(t_coarse, x_coarse, u_coarse, t_fine)

    % Dimensions
    Nd = length(t_fine); n = size(x_coarse, 1); m = size(u_coarse, 1);

    % Preallocate outputs
    x_fine = zeros(n, Nd);
    u_fine = zeros(m, Nd);

    % Interpolate each row (dimension) of x with cubic spline
    for i = 1:n
        x_fine(i, :) = interp1(t_coarse, x_coarse(i, :), t_fine, 'spline'); % Cubic interpolation
    end

    % Interpolate each row (dimension) of u with linear interpolation
    for i = 1:m
        u_fine(i, :) = interp1(t_coarse, u_coarse(i, :), t_fine, 'linear'); % Linear interpolation
    end
end

%Interpolates an input matrix
function M_fine = upsample_matrix(M_coarse, t_coarse, t_fine)
    % M_coarse: d1 × d2 × N
    % t_coarse: 1 × N
    % t_fine  : 1 × Nd (time vector corresponding to upsampling)
    
    [dim1, dim2, ~] = size(M_coarse); N_fine = length(t_fine);
    
    M_fine = zeros(dim1, dim2, N_fine);
    for i = 1:dim1
        for j = 1:dim2
            M_fine(i, j, :) = interp1(t_coarse, squeeze(M_coarse(i, j, :)), t_fine, 'linear');
        end
    end
end

function plot_xy_trajectories(trajectories, rollout_x_nom, covariance, complete_x_nom)
    % Plot all trajectories and nominal trajectory
    figure; hold on; grid on; axis equal;
    
    plot(rollout_x_nom(1, :), rollout_x_nom(2, :), 'k--', 'LineWidth', 2);

    for i = 1:length(trajectories)
        plot(trajectories{i}(1, :), trajectories{i}(2, :), 'b-', 'LineWidth', 0.5);
        plot(trajectories{i}(1, 1), trajectories{i}(2, 1), 'sg');
    end
    
    plot(complete_x_nom(1, :), complete_x_nom(2, :), 'k--', 'LineWidth', 2);
    %plot_pentagon_car(x_nom(:,end)); %final pose
    
    %plot an ellipse from which initial states are sampled
    ellipse_center = rollout_x_nom(1:2,1); % 2D center of the ellipsoid
    cov_2d = covariance(1:2, 1:2); % Extract 2D covariance
    [eig_vec, eig_val] = eig(cov_2d);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    %mu + 3*sigma will cover 99.7% of the distribution
    %std dev, sigma = sqrt(covariance)
    ellipse_boundary = 3*eig_val^(1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(ellipse_center(1) + rotated_ellipse(1, :), ...
         ellipse_center(2) + rotated_ellipse(2, :), ...
         'k-.', 'LineWidth', 1.2);

    xlabel('p_x'); ylabel('p_y');
    title('Monte Carlo Rollout Trajectories');
    legend('Nominal Trajectory','Sample Trajectories','Sample Initial states','Location','best');
    %hold off;
end

function plot_state_trajectories(trajectories, rollout_time_instances, rollout_x_nom, stateDims)
    
    if nargin < 4
        stateDims = [1 2 3];
    end

    % Plot all trajectories
    figure; %hold on; grid on; axis equal;

    for j = 1:length(stateDims)
        subplot(length(stateDims), 1, j);

        hold on; grid on;

        for i = 1:length(trajectories)
            plot(rollout_time_instances, trajectories{i}(stateDims(j), :), 'b-', 'LineWidth', 0.5);
        end
    
        plot(rollout_time_instances, rollout_x_nom(stateDims(j), :), 'k--', 'LineWidth', 2);

        xlabel('time [s]');
        ylabel(['x_{', num2str(stateDims(j)), '}'])
    end

    sgtitle('Monte Carlo Rollout Trajectories');
end

function plot_xyz_trajectories(trajectories, rollout_x_nom, covariance, complete_x_nom)
    % Plot all trajectories and nominal trajectory
    figure; view(3); hold on; grid on; axis equal;
    
    plot3(rollout_x_nom(1, :), rollout_x_nom(2, :), rollout_x_nom(3, :), 'k--', 'LineWidth', 2);

    for i = 1:length(trajectories)
        plot3(trajectories{i}(1, :), trajectories{i}(2, :), trajectories{i}(3, :), 'b-', 'LineWidth', 0.5);
        plot3(trajectories{i}(1, 1), trajectories{i}(2, 1), trajectories{i}(3, 1), 'sg');
    end
    
    plot3(complete_x_nom(1, :), complete_x_nom(2, :), complete_x_nom(3, :), 'k--', 'LineWidth', 2);
    
    % %plot an ellipse from which initial states are sampled
    % ellipse_center = rollout_x_nom(1:2,1); % 2D center of the ellipsoid
    % cov_2d = covariance(1:2, 1:2); % Extract 2D covariance
    % [eig_vec, eig_val] = eig(cov_2d);
    % 
    % theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    % %mu + 3*sigma will cover 99.7% of the distribution
    % %std dev, sigma = sqrt(covariance)
    % ellipse_boundary = 3*eig_val^(1/2) * [cos(theta); sin(theta)];
    % rotated_ellipse = eig_vec * ellipse_boundary;
    % 
    % plot(ellipse_center(1) + rotated_ellipse(1, :), ...
    %      ellipse_center(2) + rotated_ellipse(2, :), ...
    %      'k-.', 'LineWidth', 1.2);

    xlabel('p_x'); ylabel('p_y'); zlabel('p_z');
    title('Monte Carlo Rollout Trajectories');
    legend('Nominal Trajectory','Sample Trajectories','Sample Initial states','Location','best');
    %hold off;
end

function plot_input_profiles(input_profiles, rollout_time_instances, complete_u_nom, complete_time_instances)
    % Plot all trajectories and nominal trajectory
    figure;
    
    m = size(complete_u_nom,1);
    for i=1:m
        subplot(m,1,i)
        hold on; grid on; 
        for j = 1:length(input_profiles)
            plot(rollout_time_instances, input_profiles{j}(i, :), 'b-', 'LineWidth', 0.5);
        end
        plot(complete_time_instances, complete_u_nom(i, :), 'k--', 'LineWidth', 2);
        xlabel('time'); ylabel(sprintf('u_{%d}', i))

        if i==1
            title('Input history from Monte Carlo Rollouts');
        end
    end
end

% Plot error and cost metrics over time
function plot_error_metrics(errors, costs, time)
    
    figure;
    
    %error-norm
    subplot(2, 1, 1);
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
    subplot(2, 1, 2);
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