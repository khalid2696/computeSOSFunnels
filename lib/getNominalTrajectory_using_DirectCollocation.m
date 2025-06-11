function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
    getNominalTrajectory_using_DirectCollocation(cartpole_dynamics, x0, xf, T_max, N)

    %% Optimization parameters - cost function and state/input constraints
    timeHorizon_weight = 1;      % Weight for minimizing time
    controlEffort_weight = 0.01;   % Weight for control effort
    
    % State constraints
    cartPosition_limits = [-3.0, 3.0];     % Cart position bounds (m)
    cartVelocity_limits = [-5.0, 5.0];     % Cart velocity bounds (m/s)
    pendulum_angularVel_limits = [-15, 15]; % Angular velocity bounds (rad/s)
    
    % Input constraints  
    force_limits = [-20, 20];  % Horizontal force limits (N)
    
    %% Define Decision Variables
    T = sdpvar(1);              % Variable time horizon
    dt = T / (N-1);             % Adaptive time step
    
    % State Variables: [x, x_dot, theta, theta_dot]
    X = sdpvar(4, N);           % Cart position, velocity, pendulum angle, angular velocity
    U = sdpvar(1, N);           % Horizontal force on cart
    
    %% Constraints and Objective
    constraints = [0.5 <= T, T <= T_max];  % Time horizon bounds
    cost = timeHorizon_weight * T;          % Penalize total time
    
    % Add control effort to cost
    for k = 1:N-1
        cost = cost + controlEffort_weight * U(k)^2;
    end
    
    %% Dynamics Constraints (Collocation)
    for k = 1:N-1
        xk = X(:, k);
        uk = U(k);
        
        % Compute nonlinear dynamics at current and next points
        f_k = cartpole_dynamics(xk, uk);
        f_k_next = cartpole_dynamics(X(:, k+1), U(k));
        
        % Trapezoidal integration with variable time step
        x_next = xk + (dt/2) * (f_k + f_k_next);
        
        constraints = [constraints, X(:, k+1) == x_next];
    end
    
    %% Boundary Conditions
    % Initial and final state constraints
    constraints = [constraints, X(:, 1) == x0];
    constraints = [constraints, X(:, end) == xf];
    
    %% State Constraints
    % Cart position limits
    constraints = [constraints, cartPosition_limits(1) <= X(1, :)];
    constraints = [constraints, X(1, :) <= cartPosition_limits(2)];
    
    % Cart velocity limits  
    constraints = [constraints, cartVelocity_limits(1) <= X(2, :)];
    constraints = [constraints, X(2, :) <= cartVelocity_limits(2)];
    
    % Pendulum angle: no explicit limits (allow full rotation for swing-up)
    % But could add: 
    
    %constraints = [constraints, -2*pi <= X(3,:), X(3,:) <= 2*pi]; %for swing up and swing down
    constraints = [constraints, 0.9*pi <= X(3,:), X(3,:) <= 1.1*pi]; %for balancing at the top
    %constraints = [constraints, -0.1*pi <= X(3,:), X(3,:) <= 0.1*pi]; %for hanging down position
    
    % Angular velocity limits
    constraints = [constraints, pendulum_angularVel_limits(1) <= X(4, :)];
    constraints = [constraints, X(4, :) <= pendulum_angularVel_limits(2)];
    
    %% Control Input Constraints
    constraints = [constraints, force_limits(1) <= U];
    constraints = [constraints, U <= force_limits(2)];
    
    % Terminal control constraint (balanced upright requires zero force)
    constraints = [constraints, U(end) == 0];
    
    %% Additional Path Constraints (Optional)
    % Smooth control - penalize control rate changes
    % for k = 1:N-2
    %     constraints = [constraints, (U(k+1) - U(k))^2 <= max_control_rate^2];
    % end
    
    %% Solve using IPOPT
    options = sdpsettings('solver', 'ipopt', 'verbose', 1);
    
    % IPOPT-specific options for better convergence
    options.ipopt.max_iter = 3000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;
    
    diagnostics = optimize(constraints, cost, options);
    
    %% Extract Solutions
    if diagnostics.problem == 0
        fprintf('Cartpole trajectory optimization successful!\n');
        
        T_opt = value(T);
        x_opt = value(X);  % [x; x_dot; theta; theta_dot] 
        u_opt = value(U);  % [F] - horizontal force
        time_instances_opt = linspace(0, T_opt, N);
        cost_opt = value(cost);
        
        % Display results
        fprintf('Optimal time horizon: %.3f seconds\n', T_opt);
        %fprintf('Final cost: %.6f\n', cost_opt);
        fprintf('Cart final position: %.3f m\n', x_opt(1, end));
        fprintf('Pendulum final angle: %.3f rad (%.1f deg)\n', x_opt(3, end), rad2deg(x_opt(3, end)));
        
    else
        % Optimization failed
        fprintf('Cartpole trajectory optimization failed!\n');
        fprintf('Solver status: %s\n', diagnostics.info);
        
        T_opt = NaN;
        x_opt = NaN;
        u_opt = NaN;
        time_instances_opt = NaN;
        cost_opt = NaN;
    end
end