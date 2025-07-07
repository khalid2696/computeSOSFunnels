function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
    getNominalTrajectory_using_DirectCollocation_EnergyPumping(cartpole_dynamics, x0, xf, T_max, N)

    %% Optimization parameters - focus on smoothness only
    timeHorizon_weight = 0.1;
    controlEffort_weight = 0.001;
    controlSmooth_weight = 1.0;      % Penalize control derivatives
    stateSmooth_weight = 0.1;        % Penalize state derivatives
    
    % Tighter input constraints to force energy pumping strategy
    force_limits = [-6, 6];          % Even tighter to necessitate pumping
    
    % State constraints
    cartPosition_limits = [-2.0, 2.0];
    cartVelocity_limits = [-5.0, 5.0];
    pendulum_angularVel_limits = [-20, 20];
    
    %% Define Decision Variables
    T = sdpvar(1);
    dt = T / (N-1);
    
    X = sdpvar(4, N);    % [x, x_dot, theta, theta_dot]
    U = sdpvar(1, N);    % Horizontal force
    
    %% Cost Function (NO energy terms - pure smoothness)
    constraints = [2.0 <= T, T <= T_max];  % Allow sufficient time for pumping
    cost = timeHorizon_weight * T;
    
    % Control effort
    for k = 1:N
        cost = cost + controlEffort_weight * U(k)^2;
    end
    
    % Control smoothness
    for k = 1:N-1
        u_derivative = (U(k+1) - U(k)) / dt;
        cost = cost + controlSmooth_weight * u_derivative^2;
    end
    
    % State smoothness
    for k = 1:N-1
        state_derivative = (X(:,k+1) - X(:,k)) / dt;
        cost = cost + stateSmooth_weight * sum(state_derivative.^2);
    end
    
    %% Hermite-Simpson Collocation
    for k = 1:N-1
        xk = X(:, k);
        xk1 = X(:, k+1);
        uk = U(k);
        uk1 = U(k+1);
        
        % Midpoint state and control
        x_mid = 0.5 * (xk + xk1);
        u_mid = 0.5 * (uk + uk1);
        
        % Dynamics at endpoints and midpoint
        f_k = cartpole_dynamics(xk, uk);
        f_k1 = cartpole_dynamics(xk1, uk1);
        f_mid = cartpole_dynamics(x_mid, u_mid);
        
        % Hermite-Simpson integration
        x_next = xk + (dt/6) * (f_k + 4*f_mid + f_k1);
        
        constraints = [constraints, xk1 == x_next];
    end
    
    %% Boundary Conditions
    constraints = [constraints, X(:, 1) == x0];
    constraints = [constraints, X(:, end) == xf];
    
    %% ENERGY PUMPING CONSTRAINTS
    % These constraints enforce energy pumping behavior based on physics
    
    % Energy pumping strategy: Force should be applied to accelerate cart
    % in the direction that increases pendulum energy
    
    for k = 1:round(N*0.8)  % Apply pumping constraints for first 80% of trajectory
        theta_k = X(3, k);
        theta_dot_k = X(4, k);
        x_dot_k = X(2, k);
        u_k = U(k);
        
        % Energy pumping logic based on pendulum phase
        % When pendulum swings toward one side, cart should accelerate to add energy
        
        % Method 1: Velocity-based energy pumping
        % Force should have same sign as (sin(theta) * theta_dot)
        % This maximizes power transfer: P = F * v_cart_effective
        energy_pumping_signal = sin(theta_k) * theta_dot_k;
        
        % Create pumping constraints with slack variables for robustness
        epsilon = 0.1;  % Tolerance for constraint violation
        
        % If energy pumping signal is positive, force should be positive (within limits)
        % If energy pumping signal is negative, force should be negative (within limits)
        % Use complementarity-style constraints with binary logic approximation
        
        % Smooth approximation of sign function using tanh
        pumping_direction = tanh(5 * energy_pumping_signal);  % Smooth sign function
        
        % Force should align with pumping direction when energy is low
        current_energy = 0.5 * theta_dot_k^2 - cos(theta_k) + 2;  % +2 to make positive
        energy_threshold = 1.5;  % Below this, apply pumping
        
        % Only apply pumping when energy is below threshold
        energy_deficit = max(0, energy_threshold - current_energy);
        pumping_strength = energy_deficit / energy_threshold;  % 0 to 1
        
        % Energy pumping constraint: Force should align with pumping direction
        % when energy is low and angular velocity is significant
        angular_vel_threshold = 1.0;  % rad/s
        is_moving = tanh(abs(theta_dot_k) / angular_vel_threshold);
        
        % Combined pumping constraint
        desired_force_direction = pumping_strength * is_moving * pumping_direction;
        
        % Soft constraint: force should be in the energy pumping direction
        force_alignment = u_k * desired_force_direction;
        
        % Only enforce when pumping is beneficial (significant motion and low energy)
        min_alignment = -0.5 * pumping_strength * is_moving;
        constraints = [constraints, force_alignment >= min_alignment];
    end
    
    %% Alternative Energy Pumping Constraints (Method 2)
    % Enforce energy increase during specific phases
    energy_increase_phases = round(linspace(2, round(N*0.6), 5));  % Check at 5 points
    
    for i = 1:length(energy_increase_phases)-1
        k1 = energy_increase_phases(i);
        k2 = energy_increase_phases(i+1);
        
        % Calculate energy at both points
        E1 = 0.5 * X(4,k1)^2 - cos(X(3,k1));
        E2 = 0.5 * X(4,k2)^2 - cos(X(3,k2));
        
        % Energy should increase (or not decrease too much)
        min_energy_increase = -0.2;  % Allow small decreases
        constraints = [constraints, E2 >= E1 + min_energy_increase];
    end
    
    %% Resonant Frequency Pumping Constraints (Method 3)
    % Force pendulum to oscillate near its natural frequency during pumping phase
    natural_freq = sqrt(9.81 / 1);  % Assuming L=1m, g=9.81
    
    % During pumping phase, encourage oscillations near natural frequency
    pumping_phase_end = round(N * 0.7);
    for k = 2:pumping_phase_end-1
        theta_k = X(3, k);
        theta_prev = X(3, k-1);
        theta_next = X(3, k+1);
        
        % Second derivative of angle (angular acceleration)
        theta_ddot = (theta_next - 2*theta_k + theta_prev) / dt^2;
        
        % For small angles, natural oscillation: theta_ddot ≈ -w²*theta
        % But for swing-up, we want large amplitude oscillations
        % Constraint: angular acceleration should not be too large
        max_angular_accel = 50;  % rad/s²
        constraints = [constraints, -max_angular_accel <= theta_ddot <= max_angular_accel];
    end
    
    %% Standard State Constraints
    % Cart position limits
    constraints = [constraints, cartPosition_limits(1) <= X(1, :)];
    constraints = [constraints, X(1, :) <= cartPosition_limits(2)];
    
    % Cart velocity limits  
    constraints = [constraints, cartVelocity_limits(1) <= X(2, :)];
    constraints = [constraints, X(2, :) <= cartVelocity_limits(2)];
    
    % Allow full pendulum rotation - NO angle constraints for swing-up
    
    % Angular velocity limits
    constraints = [constraints, pendulum_angularVel_limits(1) <= X(4, :)];
    constraints = [constraints, X(4, :) <= pendulum_angularVel_limits(2)];
    
    %% Control Input Constraints
    constraints = [constraints, force_limits(1) <= U];
    constraints = [constraints, U <= force_limits(2)];
    
    % Control rate limits for smoothness
    max_control_rate = 30;  % N/s
    for k = 1:N-1
        du_dt = (U(k+1) - U(k)) / dt;
        constraints = [constraints, -max_control_rate <= du_dt <= max_control_rate];
    end
    
    % Terminal constraints - small force allowed for final balancing
    constraints = [constraints, -1.0 <= U(end) <= 1.0];
    
    %% Phase-Specific Constraints
    % Divide trajectory into phases: pumping phase and balancing phase
    pumping_phase = 1:round(N*0.7);
    balancing_phase = round(N*0.7)+1:N;
    
    % In balancing phase, pendulum should be moving toward upright
    for k = balancing_phase
        if k < N
            theta_k = X(3, k);
            theta_next = X(3, k+1);
            
            % Angle should be approaching π (upright)
            angle_error_k = abs(theta_k - pi);
            angle_error_next = abs(theta_next - pi);
            
            % Soft constraint: angle error should not increase dramatically
            constraints = [constraints, angle_error_next <= angle_error_k + 0.1];
        end
    end
    
    %% Solve with IPOPT
    options = sdpsettings('solver', 'ipopt', 'verbose', 1);
    
    options.ipopt.max_iter = 5000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.nlp_scaling_method = 'gradient-based';
    
    diagnostics = optimize(constraints, cost, options);
    
    %% Extract and Analyze Results
    if diagnostics.problem == 0
        fprintf('Energy pumping trajectory optimization successful!\n');
        
        T_opt = value(T);
        x_opt = value(X);
        u_opt = value(U);
        time_instances_opt = linspace(0, T_opt, N);
        cost_opt = value(cost);
        
        % Analyze energy pumping effectiveness
        analyzeEnergyPumpingConstraints(time_instances_opt, x_opt, u_opt);
        
        fprintf('Optimal time horizon: %.3f seconds\n', T_opt);
        fprintf('Cart final position: %.3f m\n', x_opt(1, end));
        fprintf('Pendulum final angle: %.3f rad (%.1f deg)\n', x_opt(3, end), rad2deg(x_opt(3, end)));
        
    else
        fprintf('Energy pumping trajectory optimization failed!\n');
        fprintf('Solver status: %s\n', diagnostics.info);
        
        T_opt = NaN; x_opt = NaN; u_opt = NaN;
        time_instances_opt = NaN; cost_opt = NaN;
    end
end

%% Analysis function for energy pumping constraints
function analyzeEnergyPumpingConstraints(time_instances, x_opt, u_opt)
    
    figure('Position', [100, 100, 1400, 800]);
    
    % Extract states
    x_pos = x_opt(1,:);
    x_vel = x_opt(2,:);
    theta = x_opt(3,:);
    theta_dot = x_opt(4,:);
    
    % Calculate energy components
    kinetic_energy = 0.5 * (x_vel.^2 + theta_dot.^2);
    potential_energy = -cos(theta);
    total_energy = kinetic_energy + potential_energy;
    
    % Calculate energy pumping signal
    energy_pumping_signal = sin(theta) .* theta_dot;
    
    % Power analysis: P = F * v_effective
    power_input = u_opt .* (x_vel + sin(theta) .* theta_dot);  % Simplified power
    
    subplot(2,3,1);
    plot(time_instances, total_energy, 'b-', 'LineWidth', 2); hold on;
    plot(time_instances, kinetic_energy, 'r--', 'LineWidth', 1);
    plot(time_instances, potential_energy, 'g--', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Energy');
    title('Energy Evolution with Pumping Constraints');
    legend('Total', 'Kinetic', 'Potential', 'Location', 'best');
    grid on;
    
    subplot(2,3,2);
    plot(time_instances, u_opt, 'k-', 'LineWidth', 2); hold on;
    plot(time_instances, 2*energy_pumping_signal, 'r--', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Force (N) / Signal');
    title('Control vs Energy Pumping Signal');
    legend('Control Force', '2×Pumping Signal', 'Location', 'best');
    grid on;
    
    subplot(2,3,3);
    plot(time_instances, rad2deg(theta), 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angle (deg)');
    title('Pendulum Angle');
    yline(180, 'r--', 'Target');
    grid on;
    
    subplot(2,3,4);
    plot(time_instances, power_input, 'g-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Power (W)');
    title('Power Input to System');
    grid on;
    
    subplot(2,3,5);
    % Phase portrait
    plot(rad2deg(theta), theta_dot, 'b-', 'LineWidth', 2); hold on;
    plot(rad2deg(theta(1)), theta_dot(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(rad2deg(theta(end)), theta_dot(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('Angle (deg)'); ylabel('Angular Velocity (rad/s)');
    title('Phase Portrait');
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
    grid on;
    
    subplot(2,3,6);
    % Force-energy correlation
    scatter(energy_pumping_signal, u_opt, 50, time_instances, 'filled');
    xlabel('Energy Pumping Signal'); ylabel('Control Force (N)');
    title('Force vs Pumping Signal (colored by time)');
    colorbar; colormap('viridis');
    grid on;
    
    sgtitle('Energy Pumping Constraint Analysis');
    
    % Print pumping effectiveness metrics
    energy_increase = total_energy(end) - total_energy(1);
    max_energy = max(total_energy);
    fprintf('\nEnergy Pumping Analysis:\n');
    fprintf('Initial energy: %.3f\n', total_energy(1));
    fprintf('Final energy: %.3f\n', total_energy(end));
    fprintf('Maximum energy reached: %.3f\n', max_energy);
    fprintf('Energy increase: %.3f\n', energy_increase);
    
    % Correlation between force and pumping signal
    correlation = corrcoef(energy_pumping_signal, u_opt);
    fprintf('Force-pumping signal correlation: %.3f\n', correlation(1,2));
end