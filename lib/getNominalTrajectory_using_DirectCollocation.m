function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
getNominalTrajectory_using_DirectCollocation(dynamics, x0, xf, input_saturations, T_max, N)

    % Optimization variables
    x = sdpvar(3, N);    % States [x; y; theta]
    u = sdpvar(2, N);    % Control inputs [v; omega]
    T = sdpvar(1, 1);    % Time horizon

    % Time step as a function of T
    dt = T / (N - 1);

    % Cost function
    Q = diag([1, 1, 0.1]); % State cost
    R = diag([0.1, 0.1]);  % Input cost
    Q_terminal = diag([10, 10, 1]); % Terminal cost
    
    cost = T; % Minimize time horizon
    for k = 1:N-1
        %cost = cost + x(:, k)' * Q * x(:, k) + u(:, k)' * R * u(:, k);
        cost = cost + u(:, k)' * R * u(:, k); %only input cost
    end
    cost = cost + (x(:, end) - xf)'*Q_terminal*(x(:, end) - xf) + u(:, end)'*R*u(:, end);

    % Constraints
    constraints = [];

    % Initial and final state constraints
    constraints = [constraints, x(:, 1) == x0];
    constraints = [constraints, x(:, end) == xf];
    %constraints = [constraints, norm(x(:, end) - xf) <= 0.5]; %a relaxed terminal constraint
    constraints = [constraints, x(3, :) >= 0]; %constraining the orientation to be always facing "forwards" (globally)

    % Time horizon constraints
    constraints = [constraints, 0 <= T, T <= T_max];

    % Control input constraints
    constraints = [constraints, -input_saturations(1) <= u(1, :), u(1, :) <= input_saturations(1)];   % Velocity limits
    constraints = [constraints, -input_saturations(2) <= u(2, :), u(2, :) <= input_saturations(2)]; % Angular velocity limits
    
    %terminal control-input constraints
    %constraints = [constraints, u(:, 1) == zeros(size(u(:, 1)))];
    %constraints = [constraints, u(:, end) == zeros(size(u(:, end)))];

    % Dynamics constraints (collocation)
    % for k = 1:N-1
    %     f_k = dynamics(x(:, k), u(:, k));
    %     f_k1 = dynamics(x(:, k+1), u(:, k+1));
    %     x_mid = 0.5 * (x(:, k) + x(:, k+1)) + 0.5 * dt * (f_k - f_k1);
    %     constraints = [constraints, x_mid == x(:, k) + dt * f_k / 2];
    % end

    % Dynamics constraints (collocation)
    for k = 1:N-1
        f1 = dynamics(x(:, k), u(:, k));
        f2 = dynamics(x(:, k+1), u(:, k+1));
        x_next = x(:, k) + 0.5 * dt * (f1 + f2); % Trapezoidal integration
        constraints = [constraints, x(:, k+1) == x_next];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'ipopt', 'verbose', 0);
    diagnostics = optimize(constraints, cost, options);

    % Results
    if diagnostics.problem == 0
        %disp('Optimization successful!');
        x_opt = value(x);
        u_opt = value(u);
        T_opt = value(T);
        cost_opt = value(cost);
        time_instances_opt = linspace(0, T_opt, N);
        %samplingTime_opt = T_opt/(N-1); %sampling-time

    else %if the optimisation fails, assign NaN values for completeness    
        x_opt = NaN;
        u_opt = NaN;
        time_instances_opt = NaN;
        cost_opt = NaN;
        %samplingTime_opt = NaN;
        %disp('Optimization failed!');
        %disp(diagnostics.info);
    end
end