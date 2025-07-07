function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
    getNominalTrajectory_using_DirectCollocation_Smooth(cartpole_dynamics, x0, xf, T_max, N)

    %% Enhanced Optimization parameters for smooth energy pumping
    timeHorizon_weight = 0.1;        % Reduced weight - allow longer time for energy pumping
    controlEffort_weight = 0.001;    % Reduced for tighter constraints
    controlSmooth_weight = 1.0;      % NEW: Penalize control derivatives
    stateSmooth_weight = 0.1;        % NEW: Penalize state derivatives  
    energyProgress_weight = 0.5;     % NEW: Reward energy increase during swing-up
    
    % Tighter input constraints for energy pumping
    force_limits = [-8, 8];          % Reduced from [-20, 20] to force energy pumping
    
    % Relaxed state constraints to allow energy pumping oscillations
    cartPosition_limits = [-2.5, 2.5];
    cartVelocity_limits = [-6.0, 6.0];
    pendulum_angularVel_limits = [-20, 20];  % Increased for energy pumping
    
    %% Define Decision Variables
    T = sdpvar(1);
    dt = T / (N-1);
    
    X = sdpvar(4, N);    % [x, x_dot, theta, theta_dot]
    U = sdpvar(1, N);    % Horizontal force
    
    %% Enhanced Cost Function for Smoothness
    constraints = [1.0 <= T, T <= T_max];  % Allow longer time
    cost = timeHorizon_weight * T;
    
    % Control effort
    for k = 1:N
        cost = cost + controlEffort_weight * U(k)^2;
    end
    
    % % Control smoothness (penalize control derivatives)
    % for k = 1:N-1
    %     u_derivative = (U(k+1) - U(k)) / dt;
    %     cost = cost + controlSmooth_weight * u_derivative^2;
    % end
    % 
    % % State smoothness (penalize excessive state derivatives) 
    % for k = 1:N-1
    %     state_derivative = (X(:,k+1) - X(:,k)) / dt;
    %     cost = cost + stateSmooth_weight * sum(state_derivative.^2);
    % end
    % 
    % % Energy progress reward for swing-up phase
    % % Reward increasing total energy during initial swing phases
    % for k = 1:round(N*0.8)  % Only for first 80% of trajectory
    %     % Total energy = kinetic + potential
    %     theta_k = X(3,k);
    %     theta_dot_k = X(4,k);
    %     x_dot_k = X(2,k);
    % 
    %     % Simplified energy terms (pendulum + cart kinetic + pendulum potential)
    %     energy_k = 0.5 * x_dot_k^2 + 0.5 * theta_dot_k^2 - cos(theta_k);
    %     cost = cost - energyProgress_weight * energy_k;  % Negative = reward
    % end
    
    %% Hermite-Simpson Collocation (Higher Order)
    % More accurate than trapezoidal for smooth trajectories
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
    
    %% Relaxed State Constraints for Energy Pumping
    % Cart position limits
    constraints = [constraints, cartPosition_limits(1) <= X(1, :)];
    constraints = [constraints, X(1, :) <= cartPosition_limits(2)];
    
    % Cart velocity limits  
    constraints = [constraints, cartVelocity_limits(1) <= X(2, :)];
    constraints = [constraints, X(2, :) <= cartVelocity_limits(2)];
    
    % CRITICAL: Allow full pendulum rotation for energy pumping swing-up
    % Remove hard angle constraints - let the pendulum swing freely
    % constraints = [constraints, -4*pi <= X(3,:), X(3,:) <= 4*pi]; % Optional: prevent excessive wrapping
    
    % Angular velocity limits
    constraints = [constraints, pendulum_angularVel_limits(1) <= X(4, :)];
    constraints = [constraints, X(4, :) <= pendulum_angularVel_limits(2)];
    
    %% Control Input Constraints
    constraints = [constraints, force_limits(1) <= U];
    constraints = [constraints, U <= force_limits(2)];
    
    % Softer terminal control constraint
    constraints = [constraints, -0.5 <= U(end), U(end) <= 0.5];  % Small final force allowed
    
    %% Path Constraints for Smoothness
    % Limit control rate of change
    max_control_rate = 50;  % N/s - adjust based on actuator capabilities
    for k = 1:N-1
        du_dt = (U(k+1) - U(k)) / dt;
        constraints = [constraints, -max_control_rate <= du_dt <= max_control_rate];
    end
    
    %% Multiple Shooting for Better Convergence
    % Add intermediate boundary conditions to prevent local minima
    % Energy should generally increase in first half
    for k = 2:round(N/3)
        theta_k = X(3,k);
        theta_dot_k = X(4,k);
        energy_k = 0.5 * theta_dot_k^2 - cos(theta_k);
        
        theta_prev = X(3,k-1);
        theta_dot_prev = X(4,k-1);
        energy_prev = 0.5 * theta_dot_prev^2 - cos(theta_prev);
        
        % Soft constraint: energy should not decrease too much
        constraints = [constraints, energy_k >= energy_prev - 0.5];
    end
    
    %% Solve with Enhanced IPOPT Settings
    options = sdpsettings('solver', 'ipopt', 'verbose', 1);
    
    % IPOPT options for smooth trajectory optimization
    options.ipopt.max_iter = 5000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.adaptive_mu_globalization = 'obj-constr-filter';
    options.ipopt.nlp_scaling_method = 'gradient-based';
    %options.ipopt.alpha_for_y = 'primal-and-dual';
    
    % Warm start if possible
    % options.ipopt.warm_start_init_point = 'yes';
    
    diagnostics = optimize(constraints, cost, options);
    
    %% Extract Solutions
    if diagnostics.problem == 0
        fprintf('Smooth cartpole trajectory optimization successful!\n');
        
        T_opt = value(T);
        x_opt = value(X);
        u_opt = value(U);
        time_instances_opt = linspace(0, T_opt, N);
        cost_opt = value(cost);
        
        % Analyze trajectory smoothness
        control_derivatives = diff(u_opt) ./ diff(time_instances_opt);
        max_control_rate_actual = max(abs(control_derivatives));
        
        fprintf('Optimal time horizon: %.3f seconds\n', T_opt);
        fprintf('Cart final position: %.3f m\n', x_opt(1, end));
        fprintf('Pendulum final angle: %.3f rad (%.1f deg)\n', x_opt(3, end), rad2deg(x_opt(3, end)));
        fprintf('Max control rate: %.2f N/s\n', max_control_rate_actual);
        
        % Check energy pumping effectiveness
        initial_energy = 0.5 * x_opt(4,1)^2 - cos(x_opt(3,1));
        max_energy = max(0.5 * x_opt(4,:).^2 - cos(x_opt(3,:)));
        fprintf('Energy pumping ratio: %.2f\n', max_energy / abs(initial_energy));
        
    else
        fprintf('Smooth cartpole trajectory optimization failed!\n');
        fprintf('Solver status: %s\n', diagnostics.info);
        
        T_opt = NaN; x_opt = NaN; u_opt = NaN;
        time_instances_opt = NaN; cost_opt = NaN;
    end
end

% %% Additional helper function for energy analysis
% function analyzeEnergyPumping(time_instances, x_opt, u_opt)
%     % Analyze the energy pumping characteristics
% 
%     figure('Position', [100, 100, 1200, 400]);
% 
%     % Calculate total energy over time
%     theta = x_opt(3,:);
%     theta_dot = x_opt(4,:);
%     x_dot = x_opt(2,:);
% 
%     % Total energy components
%     kinetic_energy = 0.5 * (x_dot.^2 + theta_dot.^2);
%     potential_energy = -cos(theta);  % Simplified potential energy
%     total_energy = kinetic_energy + potential_energy;
% 
%     subplot(1,3,1);
%     plot(time_instances, total_energy, 'b-', 'LineWidth', 2); hold on;
%     plot(time_instances, kinetic_energy, 'r--', 'LineWidth', 1.5);
%     plot(time_instances, potential_energy, 'g--', 'LineWidth', 1.5);
%     xlabel('Time (s)'); ylabel('Energy');
%     title('Energy Evolution During Swing-Up');
%     legend('Total', 'Kinetic', 'Potential', 'Location', 'best');
%     grid on;
% 
%     subplot(1,3,2);
%     plot(time_instances, u_opt, 'k-', 'LineWidth', 2);
%     xlabel('Time (s)'); ylabel('Control Force (N)');
%     title('Control Input for Energy Pumping');
%     grid on;
% 
%     subplot(1,3,3);
%     plot(time_instances, rad2deg(theta), 'b-', 'LineWidth', 2); hold on;
%     yline(180, 'r--', 'Target Upright');
%     xlabel('Time (s)'); ylabel('Pendulum Angle (deg)');
%     title('Pendulum Angle During Swing-Up');
%     grid on;
% 
%     sgtitle('Energy Pumping Analysis');
% end