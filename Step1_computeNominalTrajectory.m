%clc; clearvars; close all
%yalmip('clear');

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
    x0 = [0; 0; 0; 0;];   %Initial state: pendulum hanging down, cart at origin
    xf = [1; 0; pi; 0;];  %Final state: pendulum balanced upright, cart at specific position
end

if exist('numTimeSteps','var')
    N = numTimeSteps;
else
    N = 50; % default number of time steps
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
%[x_nom, u_nom, time_instances, nom_trajCost, diagnostics] = getNominalTrajectory_using_DirectCollocation_Smooth(dynamicsFnHandle, x0, xf, T_max, N);
%[x_nom, u_nom, time_instances, nom_trajCost, diagnostics] = getNominalTrajectory_using_DirectCollocation_EnergyPumping(dynamicsFnHandle, x0, xf, T_max, N);
toc   

disp('Finshed computing nominal trajectory and nominal (feedforward) input tape');
disp(' ');

if diagnostics.problem ~= 0 %if the optimisation fails, relax the time-horizon and input constraints
    disp(diagnostics.info);
    error('Failed to compute a feasible nominal trajectory! Cannot proceed.. Exitting!');
end

%% save the nominal trajectory and feedforward input
save('./precomputedData/nominalTrajectory.mat', 'time_instances', 'x_nom', 'u_nom', 'dynamicsFnHandle', 'cartPoleParameters');

disp('Saved the nominal trajectory and nominal inputs to a file!');
disp(' ');

%% Additionally, plot if enabled

if drawFlag 

    animate_cartpole_trajectory(time_instances, x_nom, cartPoleParameters);

    plotPolePosition(x_nom, cartPoleParameters);
    title('Cartpole Trajectory');

    plotStateInputTrajectories(time_instances, x_nom, u_nom);
    sgtitle(sprintf('Cartpole Swing-Up Trajectory (T_f = %.2f s)', time_instances(end)));
    %sgtitle: Title for a grid of subplots
end

%clearvars; %cleanup the workspace after saving relevant data and plotting

%% Local Function definitions

function plotStateInputTrajectories(t_opt, x_opt, u_opt)
    figure('Position', [100, 100, 1200, 800]);
        
    subplot(2,3,1);
    plot(t_opt, x_opt(1,:), 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Cart Position (m)');
    title('Cart Position vs Time'); grid on;
    
    subplot(2,3,2);
    plot(t_opt, rad2deg(x_opt(3,:)), 'g-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Pendulum Angle (deg)');
    title('Pendulum Angle vs Time'); grid on;
    yline(180, 'k--', 'Upright');
    
    subplot(2,3,3);
    plot(t_opt, u_opt, 'k-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Control Force (N)');
    title('Control Input vs Time'); grid on;

    subplot(2,3,4);
    plot(t_opt, x_opt(2,:), 'r-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Cart Velocity (m/s)');
    title('Cart Velocity vs Time'); grid on;

    subplot(2,3,5);
    plot(t_opt, x_opt(4,:), 'm-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Pendulum Angular Velocity vs Time'); grid on;
    
    subplot(2,3,6);
    % Phase portrait: pendulum angle vs angular velocity
    plot(rad2deg(x_opt(3,:)), x_opt(4,:), 'b-', 'LineWidth', 2);
    xlabel('Pendulum Angle (deg)'); ylabel('Angular Velocity (rad/s)');
    title('Pendulum Phase Portrait'); grid on;
    hold on;
    plot(rad2deg(x_opt(3,1)), x_opt(4,1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(rad2deg(x_opt(3,end)), x_opt(4,end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    legend('Trajectory', 'Start', 'End');

end

% 2D Cartpole Visualization - Cart and Pole Tip Trajectory
function plotPolePosition(x_opt, params)

    figure; grid on; axis equal; hold on;
    L = params.L;  % Pendulum length
    
    % Calculate pole tip positions
    cart_x = x_opt(1,:);  % Cart x positions
    cart_y = zeros(size(cart_x));  % Cart y position (always 0)
    
    % Pole tip coordinates (pendulum length L, angle theta from vertical down)
    pole_tip_x = cart_x + L * sin(pi - x_opt(3,:));  % x = cart_x + L*sin(theta)
    pole_tip_y = L * cos(pi - x_opt(3,:));           % y = L*cos(theta) (positive up)
    
    % Plot pole tip trajectory
    plot(pole_tip_x, pole_tip_y, 'b--', 'LineWidth', 1.2); 
    
    % Plot cart trajectory (along ground)
    plot(cart_x, cart_y, 'k-', 'LineWidth', 3);
    
    % Mark start and end positions
    plot(pole_tip_x(1), pole_tip_y(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(pole_tip_x(end), pole_tip_y(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    %plot(cart_x(1), cart_y(1), 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    %plot(cart_x(end), cart_y(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Add some intermediate pole positions for visualization
    n_intermediate = 6;
    indices = round(linspace(1, length(cart_x), n_intermediate));
    for i = indices
        % Draw pole as line from cart to tip
        plot([cart_x(i), pole_tip_x(i)], [cart_y(i), pole_tip_y(i)], ...
             'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.5);
        % Draw cart as small rectangle
        rectangle('Position', [cart_x(i)-0.1, cart_y(i)-0.075, 0.2, 0.15], ...
                 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'k');
    end
    
    % Formatting
    xlabel('X Position (m)'); %ylabel('Y Position (m)');
    
    legend('Pole Tip Path', 'Cart Path', 'Start (Pole)', 'End (Pole)', ...
            'Location', 'best'); %'Start (Cart)', 'End (Cart)',
    
    % Set reasonable axis limits
    x_range = [min([cart_x, pole_tip_x]) - 0.2, max([cart_x, pole_tip_x]) + 0.2];
    y_range = [-L*1.2, L*1.2];
    xlim(x_range); ylim(y_range);
end

% Animated visualization function
function animate_cartpole_trajectory(t_opt, x_opt, params)
    % Animate the cartpole swing-up trajectory
    
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

end