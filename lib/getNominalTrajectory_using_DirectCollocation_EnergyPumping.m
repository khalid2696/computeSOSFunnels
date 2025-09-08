function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
    getNominalTrajectory_using_DirectCollocation_EnergyPumping_New(cartpole_dynamics, x0, xf, T_max, N)

    %% weights and limits
    timeHorizon_weight = 1;
    controlEffort_weight = 0.01;
    cartPosition_limits = [-3.0, 3.0];
    cartVelocity_limits = [-5.0, 5.0];
    pendulum_angularVel_limits = [-15, 15];
    force_limits = [-20, 20];

    %% Energy pumping params (tune these)
    k_ep = 20.0;     %larger k_ep --> more oscillations 
    %k_ep = 0;
    lambda_ep = 50.0; %larger lambda_ep --> more swinging up "adherance"
    %lambda_ep = 0;
    taper_power = 3;

    %% identify the operation mode
    initialPendulumAngle = x0(3);
    finalPendulumAngle = xf(3);

    if initialPendulumAngle == 0 && finalPendulumAngle == 0 %down-balance
        operationMode = 'downBalance';
    elseif initialPendulumAngle == pi && finalPendulumAngle == pi
        operationMode = 'topBalance';
    elseif initialPendulumAngle == 0 && finalPendulumAngle == pi
        operationMode = 'swingUp';
    elseif initialPendulumAngle == pi && finalPendulumAngle == 0
        operationMode = 'swingDown';
    end

    %% decision vars
    T = sdpvar(1);
    dt = T / (N-1);
    X = sdpvar(4, N);
    U = sdpvar(1, N);

    constraints = [0.5 <= T, T <= T_max];
    cost = timeHorizon_weight * T;

    % dynamics collocation + cost
    for k = 1:N-1
        xk = X(:, k);
        uk = U(k);
        f_k = cartpole_dynamics(xk, uk);
        f_k_next = cartpole_dynamics(X(:, k+1), U(k));
        x_next = xk + (dt/2) * (f_k + f_k_next);
        constraints = [constraints, X(:, k+1) == x_next];

        % control effort
        cost = cost + controlEffort_weight * U(k)^2;

        % energy pumping soft penalty
        tau_k = (k-1)/(N-1);
        s_k = 1 - tau_k^taper_power;
        U_ep = k_ep * X(4,k) * cos(X(3,k));
        cost = cost + lambda_ep * s_k * (U(k) - U_ep)^2;
    end

    % boundary & state/input limits
    constraints = [constraints, X(:,1) == x0, X(:,end) == xf];
    constraints = [constraints, cartPosition_limits(1) <= X(1,:), X(1,:) <= cartPosition_limits(2)];
    constraints = [constraints, cartVelocity_limits(1) <= X(2,:), X(2,:) <= cartVelocity_limits(2)];
    constraints = [constraints, pendulum_angularVel_limits(1) <= X(4,:), X(4,:) <= pendulum_angularVel_limits(2)];
    constraints = [constraints, force_limits(1) <= U, U <= force_limits(2)];
    constraints = [constraints, U(end) == 0];



    options = sdpsettings('solver', 'ipopt', 'verbose', 1);
    options.ipopt.max_iter = 3000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;

    diagnostics = optimize(constraints, cost, options);

    if diagnostics.problem == 0
        T_opt = value(T);
        x_opt = value(X);
        u_opt = value(U);
        time_instances_opt = linspace(0, T_opt, N);
        cost_opt = value(cost);
    else
        T_opt = NaN; x_opt = NaN; u_opt = NaN; time_instances_opt = NaN; cost_opt = NaN;
    end
end
