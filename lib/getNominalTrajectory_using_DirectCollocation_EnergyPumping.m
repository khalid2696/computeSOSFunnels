function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
    getNominalTrajectory_using_DirectCollocation_EnergyPumping(cartpole_dynamics, x0, xf, T_max, N)

    %% Energy Pumping Parameters
    k_pump = 5.0;                    % Energy pumping gain
    energy_threshold = 0.8;          % Switch to balancing when energy > threshold * target_energy
    balancing_region = pi/6;         % Angular region around upright for balancing (rad)
    
    % Physical parameters (assume these are available or pass them in)
    m = 0.5;    % Pendulum mass (kg) - adjust based on your cartPoleParameters
    g = 9.81;   % Gravity (m/s²)
    L = 0.5;    % Pendulum length (m) - adjust based on your cartPoleParameters
    target_energy = m * g * L;       % Energy needed to reach upright position
    
    %% Optimization parameters
    timeHorizon_weight = 1;
    controlEffort_weight = 0.01;
    energy_pumping_weight = 10;      % Weight for energy pumping objective
    
    % State constraints
    cartPosition_limits = [-3.0, 3.0];
    cartVelocity_limits = [-5.0, 5.0];
    pendulum_angularVel_limits = [-15, 15];
    
    % Input constraints  
    force_limits = [-20, 20];
    
    %% Define Decision Variables
    T = sdpvar(1);
    dt = T / (N-1);
    
    X = sdpvar(4, N);   % [x, x_dot, theta, theta_dot]
    U = sdpvar(1, N);   % Control force
    
    %% APPROACH 1: Energy Pumping as Soft Constraint in Cost Function
    constraints = [0.5 <= T, T <= T_max];
    cost = timeHorizon_weight * T;
    
    % Add control effort
    for k = 1:N-1
        cost = cost + controlEffort_weight * U(k)^2;
    end
    
    % Add energy pumping objective - simplified approach
    for k = 1:N-1
        theta = X(3, k);
        omega = X(4, k);
        
        % Desired energy pumping force
        u_desired = k_pump * omega * cos(theta);
        
        % Current pendulum energy (relative to hanging down position)
        current_energy = 0.5 * m * L^2 * omega^2 + m * g * L * (1 - cos(theta));
        
        % Simple weighting based on energy deficit (linear, no exponentials)
        energy_deficit = target_energy - current_energy;
        energy_weight = energy_deficit / target_energy;  % This will be 0-1 when energy is below target
        
        % Additional weighting based on distance from upright (quadratic function)
        angle_deviation = (theta - pi)^2;  % 0 when at upright, larger when away
        angle_weight = angle_deviation / (pi^2);  % Normalize to 0-1 range
        
        % Combined weighting (both must be significant for energy pumping to be active)
        pumping_weight = energy_weight * angle_weight;
        
        % Add energy pumping term to cost
        energy_pumping_error = (U(k) - u_desired)^2;
        cost = cost + energy_pumping_weight * pumping_weight * energy_pumping_error;
    end
    
    %% Dynamics Constraints (same as original)
    for k = 1:N-1
        xk = X(:, k);
        uk = U(k);
        
        f_k = cartpole_dynamics(xk, uk);
        f_k_next = cartpole_dynamics(X(:, k+1), U(k));
        
        x_next = xk + (dt/2) * (f_k + f_k_next);
        constraints = [constraints, X(:, k+1) == x_next];
    end
    
    %% Boundary Conditions
    constraints = [constraints, X(:, 1) == x0];
    constraints = [constraints, X(:, end) == xf];
    
    %% State Constraints
    constraints = [constraints, cartPosition_limits(1) <= X(1, :)];
    constraints = [constraints, X(1, :) <= cartPosition_limits(2)];
    constraints = [constraints, cartVelocity_limits(1) <= X(2, :)];
    constraints = [constraints, X(2, :) <= cartVelocity_limits(2)];
    
    % Allow full rotation for swing-up phase
    constraints = [constraints, pendulum_angularVel_limits(1) <= X(4, :)];
    constraints = [constraints, X(4, :) <= pendulum_angularVel_limits(2)];
    
    %% Control Input Constraints
    constraints = [constraints, force_limits(1) <= U];
    constraints = [constraints, U <= force_limits(2)];
    
    % Terminal constraint - balanced upright
    constraints = [constraints, U(end) == 0];
    constraints = [constraints, 0.9*pi <= X(3,end), X(3,end) <= 1.1*pi];
    
    %% APPROACH 2: Alternative - Direct Energy Pumping Constraints
    % Uncomment this section for a simpler approach with direct weighting
    
    % for k = 1:N-1
    %     theta = X(3, k);
    %     omega = X(4, k);
    %     
    %     % Simple energy pumping objective without logical operations
    %     u_desired = k_pump * omega * cos(theta);
    %     
    %     % Add energy pumping term to cost (always active, but weighted by energy deficit)
    %     current_energy = 0.5 * m * L^2 * omega^2 + m * g * L * (1 - cos(theta));
    %     energy_deficit = max(0, target_energy - current_energy);
    %     energy_weight = energy_deficit / target_energy;  % 0 to 1 weighting
    %     
    %     cost = cost + energy_pumping_weight * energy_weight * (U(k) - u_desired)^2;
    % end
    
    %% Solve
    options = sdpsettings('solver', 'ipopt', 'verbose', 1);
    options.ipopt.max_iter = 3000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;
    
    diagnostics = optimize(constraints, cost, options);
    
    %% Extract Solutions
    if diagnostics.problem == 0
        fprintf('Energy pumping trajectory optimization successful!\n');
        
        T_opt = value(T);
        x_opt = value(X);
        u_opt = value(U);
        time_instances_opt = linspace(0, T_opt, N);
        cost_opt = value(cost);
        
        % Analyze energy pumping behavior
        analyze_energy_pumping_trajectory(time_instances_opt, x_opt, u_opt, k_pump, m, g, L);
        
        fprintf('Optimal time horizon: %.3f seconds\n', T_opt);
        fprintf('Cart final position: %.3f m\n', x_opt(1, end));
        fprintf('Pendulum final angle: %.3f rad (%.1f deg)\n', x_opt(3, end), rad2deg(x_opt(3, end)));
        
    else
        fprintf('Energy pumping trajectory optimization failed!\n');
        fprintf('Solver status: %s\n', diagnostics.info);
        
        T_opt = NaN;
        x_opt = NaN;
        u_opt = NaN;
        time_instances_opt = NaN;
        cost_opt = NaN;
    end
end

%% Analysis Function
function analyze_energy_pumping_trajectory(t, x_opt, u_opt, k_pump, m, g, L)
    
    theta = x_opt(3, :);
    omega = x_opt(4, :);
    
    % Calculate energy at each time step
    kinetic_energy = 0.5 * m * L^2 * omega.^2;
    potential_energy = m * g * L * (1 - cos(theta));
    total_energy = kinetic_energy + potential_energy;
    target_energy = m * g * L;
    
    % Calculate desired energy pumping force
    u_desired = k_pump * omega .* cos(theta);
    
    figure('Position', [100, 100, 1400, 900]);
    
    % Plot 1: Energy evolution
    subplot(2,3,1);
    plot(t, kinetic_energy, 'b-', 'LineWidth', 2); hold on;
    plot(t, potential_energy, 'r-', 'LineWidth', 2);
    plot(t, total_energy, 'k-', 'LineWidth', 2);
    yline(target_energy, 'g--', 'Target Energy', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Energy (J)');
    title('Energy Evolution');
    legend('Kinetic', 'Potential', 'Total', 'Target', 'Location', 'best');
    grid on;
    
    % Plot 2: Control comparison
    subplot(2,3,2);
    plot(t, u_opt, 'b-', 'LineWidth', 2); hold on;
    plot(t, u_desired, 'r--', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Control Input vs Desired Energy Pumping');
    legend('Actual U', 'Desired U_{pump}', 'Location', 'best');
    grid on;
    
    % Plot 3: Pendulum angle
    subplot(2,3,3);
    plot(t, rad2deg(theta), 'g-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angle (deg)');
    title('Pendulum Angle');
    yline(180, 'k--', 'Upright');
    grid on;
    
    % Plot 4: Angular velocity
    subplot(2,3,4);
    plot(t, omega, 'm-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Pendulum Angular Velocity');
    grid on;
    
    % Plot 5: Energy pumping effectiveness
    subplot(2,3,5);
    energy_ratio = total_energy / target_energy;
    plot(t, energy_ratio * 100, 'k-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Energy Ratio (%)');
    title('Energy as % of Target');
    yline(100, 'r--', '100%');
    grid on;
    
    % Plot 6: Phase portrait
    subplot(2,3,6);
    plot(rad2deg(theta), omega, 'b-', 'LineWidth', 2);
    xlabel('Angle (deg)'); ylabel('Angular Velocity (rad/s)');
    title('Phase Portrait');
    hold on;
    plot(rad2deg(theta(1)), omega(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(rad2deg(theta(end)), omega(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    grid on;
    
    sgtitle('Energy Pumping Trajectory Analysis');
end

%% Alternative: Pure Energy Pumping with IVP Integration
function [x_opt, u_opt, time_instances_opt] = energy_pumping_IVP_approach(cartpole_dynamics, x0, xf, T_max, dt)
    % This is a simpler approach using just energy pumping control law
    % integrated with an IVP solver, then post-processed for the framework
    
    % Parameters
    k_pump = 5.0;
    m = 0.5; g = 9.81; L = 0.5;
    target_energy = m * g * L;
    energy_threshold = 0.95;
    
    % Time vector
    t = 0:dt:T_max;
    N = length(t);
    
    % Initialize
    x = zeros(4, N);
    u = zeros(1, N);
    x(:, 1) = x0;
    
    % Simulate with energy pumping control
    for k = 1:N-1
        theta = x(3, k);
        omega = x(4, k);
        
        % Current energy
        current_energy = 0.5 * m * L^2 * omega^2 + m * g * L * (1 - cos(theta));
        energy_ratio = current_energy / target_energy;
        
        % Control law
        if energy_ratio < energy_threshold && abs(theta - pi) > pi/6
            % Energy pumping phase
            u(k) = k_pump * omega * cos(theta);
        else
            % Balancing phase (simple PD control)
            kp = 50; kd = 10;
            u(k) = kp * (pi - theta) - kd * omega;
        end
        
        % Integrate dynamics
        x_dot = cartpole_dynamics(x(:, k), u(k));
        x(:, k+1) = x(:, k) + dt * x_dot;
    end
    
    % Return in same format as direct collocation
    x_opt = x;
    u_opt = u;
    time_instances_opt = t;
end