function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
getNominalTrajectory_using_DirectCollocation(x0, xf, T_max, N)
    
    %% Quadrotor parameters
    m = 0.5; % Mass (kg)
    g = 9.81; % acc due to gravity (m/s²)
    J = diag([0.01, 0.01, 0.018]); % Inertia matrix (kg-m^2)

    %% Optimisation parameters - cost function and state/input constraints
    timeHorizon_weight = 3;
    controlEffort_weight = 0.1;
    
    thrustSaturation = [0 15];
    momentSaturation = [-2 2];
    yawSaturation = [-0.25 0.25]; %heading angle corresponding counterpart in quaternion representation

    % Decision variable for time horizon
    T = sdpvar(1); 
    dt = T / (N-1); % Adaptive time step
    
    % State Variables: [pos, vel, quat, omega]
    X = sdpvar(13, N); % Position (3), velocity (3), quaternion (4), angular velocity (3)
    U = sdpvar(4, N); % Thrust (1) and torques (3)
    
    % Constraints and Objective
    constraints = [0.1 <= T, T <= T_max]; % Time horizon bounds
    cost = timeHorizon_weight * T; % Penalize total time taken
    
    %Objective cost
    for k = 1:N-1
        % Cost function: time + tracking error + control effort
        %cost = cost + norm(X(1:3, k) - pf, 2)^2 + 0.1 * norm(U(:, k), 2)^2;
        
        %cost function: time + control effort
        cost = cost + controlEffort_weight * norm(U(:, k), 2)^2;
    end
    
    % Dynamics constraints (collocation)
    for k = 1:N-1
        xk = X(:, k);
        uk = U(:, k);
        
        % Compute nonlinear dynamics
        f_k = quadrotor_dynamics(xk, uk, m, g, J);
        f_k_next = quadrotor_dynamics(X(:, k+1), uk, m, g, J);
    
        %x_next = xk + dt*f_k; %Euler integration
        x_next = xk + (dt/2) * (f_k + f_k_next); %trapezoidal integration
        
        % Trapezoidal integration with variable time step
        constraints = [constraints, X(:, k+1) == x_next];
    end
    
    % Initial and final state constraints
    constraints = [constraints, X(:,1) == x0, X(:,end) == xf];
    constraints = [constraints, yawSaturation(1) <= X(10,:), X(10,:) <= yawSaturation(2)]; % Heading angle
    
    % Control input constraints -- thrust and moments 
    % Rotor speeds can be computed using the mixer model with parameters: 
    % thrust coefficient (k), counter-moment coefficient (d) and arm length (l)
    constraints = [constraints, thrustSaturation(1) <= U(1,:), U(1,:) <= thrustSaturation(2)]; % Thrust
    constraints = [constraints, momentSaturation(1) <= U(2:4,:), U(2:4,:) <= momentSaturation(2)]; % Torques
    constraints = [constraints, U(1,end) == m*g]; %hover at terminal state
    
    % Solve using IPOPT
    options = sdpsettings('solver', 'ipopt', 'verbose', 0);
    diagnostics = optimize(constraints, cost, options);
    
    %% Extract solutions
    
    if diagnostics.problem == 0
        %disp('Optimization successful!');
        T_opt = value(T); %scalar - time
        x_opt = value(X); %x,y,z,v_x,v_y,v_z,q0,q1,q2,q3,p,q,r
        u_opt = value(U); %T, M_p, M_q, M_r
    
        time_instances_opt = linspace(0,T_opt,N); %may have to take transpose
        cost_opt = value(cost);
        %samplingTime_opt = T_opt/(N-1); %sampling-time
    
    else %if the optimisation fails, assign NaN values for completeness    
        time_instances_opt = NaN;
        x_opt = NaN;
        u_opt = NaN;
        time_instances_opt = NaN;
        cost_opt = NaN;
        %samplingTime_opt = NaN;
        %disp('Optimization failed!');
        %disp(diagnostics.info);
    end
end

function xdot = quadrotor_dynamics(x, u, m, g, J)
    % State variables
    p = x(1:3); % Position
    v = x(4:6); % Velocity
    q = x(7:10); % Quaternion
    omega = x(11:13); % Angular velocity
    
    % Control inputs
    T = u(1); % Thrust [N]
    tau = u(2:4); % Torques [N-m]
    
    % Quaternion components
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    % Define the symbolic rotation matrix R(q)
    R = [1 - 2*(q3^2 + q2^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
 

    % Translational dynamics
    v_dot = [0; 0; -g] + (R * [0; 0; T]) / m;
    
    % Rotational dynamics
    omega_dot = J \ (tau - cross(omega, J * omega));
    
    % Quaternion derivative
    q_dot = 0.5 * quatmultiply(q', [0; omega]')';
    
    % Return full state derivative
    xdot = [v; v_dot; q_dot; omega_dot];
end