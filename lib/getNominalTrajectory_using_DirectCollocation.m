function [x_opt, u_opt, time_instances_opt, cost_opt, diagnostics] = ...
getNominalTrajectory_using_DirectCollocation(quadrotor_dynamics, x0, xf, T_max, N, quadParameters)

    %% extracting quadrotor model parameters
    m = quadParameters.m; g = quadParameters.g;

    %% Optimisation parameters - cost function and state/input constraints
    timeHorizon_weight = 0.75; %1
    controlEffort_weight = 0.25;
    controlSmoothnessWeight = 20;  %10    % NEW: Penalize control derivatives
    velocitySmoothnessWeight = 20; %20    % NEW: Penalize state derivatives
    
    %state constraints
    rollSaturation = [-pi/6 pi/6]; %small roll and pitch angles
    pitchSaturation = [-pi/6 pi/6];
    yawSaturation = [-pi/6 pi/6]; %heading angle

    %input saturations
    Weight = m*g;
    thrustSaturation = [0 2*Weight]; %considering a T/W ratio of 2
    momentSaturation = [-1 1];
    
    %% Define Decision variables 
    T = sdpvar(1); %for time horizon
    %T = T_max;
    dt = T / (N-1); % Adaptive time step

    % State Variables: [pos, vel, attitude, omega]
    X = sdpvar(12, N); % Position (3), velocity (3), Euler angles (3), angular velocity (3)
    U = sdpvar(4, N); % Thrust (1) and torques (3)  
    
    %% Cost Objective
    % 
    % cost = timeHorizon_weight * T; % Penalize total time taken
    % 
    % for k = 1:N-1
    %     % Cost function: time + tracking error + control effort
    %     %cost = cost + norm(X(1:3, k) - pf, 2)^2 + 0.1 * norm(U(:, k), 2)^2;
    % 
    %     %cost function: time + control effort
    %     cost = cost + controlEffort_weight * norm(U(:, k), 2)^2;
    % end

    %% Cost Objective with smoothness regularization
    %cost = 0;
    cost = timeHorizon_weight * T; % Base time penalty
    
    for k = 1:N-1
        % Basic control effort penalty
        cost = cost + controlEffort_weight * (U(1,k)^2 + U(2,k)^2 + U(3,k)^2 + U(4,k)^2);
        
        % NEW: Control smoothness penalty - penalize rapid changes in inputs
        if k > 1
            thrust_change = U(1,k) - U(1,k-1);
            cost = cost + controlSmoothnessWeight * thrust_change^2;
            
            % Penalize changes in moments
            for i = 2:4
                moment_change = U(i,k) - U(i,k-1);
                cost = cost + controlSmoothnessWeight * moment_change^2;
            end
        end
        
        % NEW: State smoothness penalty - penalize rapid changes in velocities and angular rates
        if k > 1
            for i = 4:6  % linear velocity components
                vel_change = X(i,k) - X(i,k-1);
                cost = cost + velocitySmoothnessWeight * vel_change^2;
            end

            for i = 10:12  % angular velocity components
                vel_change = X(i,k) - X(i,k-1);
                cost = cost + velocitySmoothnessWeight * vel_change^2;
            end
        end
    end
    
    %% Constraints

    constraints = [];
    constraints = [0.1 <= T, T <= T_max]; % Time horizon bounds

    % Dynamics constraints (collocation)
    for k = 1:N-1
        xk = X(:, k);
        uk = U(:, k);
        
        % Compute nonlinear dynamics

        %Euler integration
        %f_k = quadrotor_dynamics(xk, uk);
        %x_next = xk + dt*f_k;
        
        %Trapezoidal integration
        f_k = quadrotor_dynamics(xk, uk);
        f_k_next = quadrotor_dynamics(X(:, k+1), uk);
        x_next = xk + (dt/2) * (f_k + f_k_next);

        %RK4 integration
        %k1 = quadrotor_dynamics(xk, uk);
        %k2 = quadrotor_dynamics(xk + 0.5 * dt * k1, uk);
        %k3 = quadrotor_dynamics(xk + 0.5 * dt * k2, uk);
        %k4 = quadrotor_dynamics(xk + dt * k3, uk);
        %x_next = xk + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4); % Update state
        
        % Enforce dynamics constraints
        constraints = [constraints, X(:, k+1) == x_next];
    end

    % Initial and final state constraints
    constraints = [constraints, X(:,1) == x0, X(:,end) == xf];
    
    %state constraints
    % height should be greater than zero
    constraints = [constraints, 0 <= X(3,:)]; 
    
    %small roll and pitch angles
    constraints = [constraints, rollSaturation(1) <= X(7,:), X(7,:) <= rollSaturation(2)];
    constraints = [constraints, pitchSaturation(1) <= X(8,:), X(8,:) <= pitchSaturation(2)];
    
    %heading angle (yaw) constrained to not deviate too much from 'North'
    constraints = [constraints, yawSaturation(1) <= X(9,:), X(9,:) <= yawSaturation(2)];
    
    % Control input constraints -- thrust and moments 
    % Rotor speeds can be computed using the mixer model with parameters: 
    % thrust coefficient (k), counter-moment coefficient (d) and arm length (l)
    constraints = [constraints, thrustSaturation(1) <= U(1,:), U(1,:) <= thrustSaturation(2)]; % Thrust
    constraints = [constraints, momentSaturation(1) <= U(2,:), U(2,:) <= momentSaturation(2)]; % Torques
    constraints = [constraints, momentSaturation(1) <= U(3,:), U(3,:) <= momentSaturation(2)]; % Torques
    constraints = [constraints, momentSaturation(1) <= U(4,:), U(4,:) <= momentSaturation(2)]; % Torques

    %hover at terminal state
    constraints = [constraints, U(1,end:end) == Weight]; %hover at terminal state
    constraints = [constraints, U(2,end:end) == 0];
    constraints = [constraints, U(3,end:end) == 0];
    constraints = [constraints, U(4,end:end) == 0];
    
    % %constraints on jerky inputs
    % for k=2:N
    %     thrust_change = U(1,k) - U(1,k-1);
    %     constraints = [constraints, -2 <= thrust_change/dt, thrust_change/dt <= 2]; 
    % 
    %         % Penalize changes in moments
    %         for i = 2:4
    %             moment_change = U(i,k) - U(i,k-1);
    %             constraints = [constraints, -2 <= moment_change/dt, moment_change/dt <= 2];
    %         end
    % end
    
    %% Solve using IPOPT
    %options = sdpsettings('solver', 'ipopt', 'verbose', 2);
    options = sdpsettings('solver', 'ipopt', 'verbose', 0);
    options.ipopt.max_iter = 5000; % Set the desired tolerance
    options.ipopt.tol = 1e-6; % Set the desired tolerance
    options.ipopt.acceptable_tol = 1e-4; % Set the acceptable tolerance 
    options.ipopt.acceptable_iter = 100; % Number of iterations for acceptable termination
    
    %numVars = getvariables(constraints);
    %length(unique(numVars))
    
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

% % Quadrotor dynamics
% function f = quadrotor_dynamics(x, u, quadParameters)
%     % Numerical version with specific parameter values
% 
%     %extracting quadrotor model parameters
%     m = quadParameters.m; g = quadParameters.g;
%     Jxx = quadParameters.J(1); Jyy = quadParameters.J(2); Jzz = quadParameters.J(3);
% 
%     %assigning state variables for ease of usage
%     % State: x = [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
% 
%     %px = x(1); py = x(2); pz = x(3);
%     vx = x(4); vy = x(5); vz = x(6);
%     phi = x(7); theta = x(8); psi = x(9);
%     p = x(10); q = x(11); r = x(12);
% 
%     %assigning input variables for ease of usage
%     % Input: u = [T; Mx; My; Mz]
%     T = u(1); Mp = u(2); Mq = u(3); Mr = u(4);
% 
%     % Trigonometric shortcuts
%     c_phi = cos(phi); s_phi = sin(phi);
%     c_theta = cos(theta); s_theta = sin(theta);
%     c_psi = cos(psi); s_psi = sin(psi);
%     t_theta = tan(theta);
%     sec_theta = sec(theta);
% 
%     % Quadrotor dynamics 
%     % Assumptions: no aerodynamic drag and gyroscopic coupling due to rotor inertia)
%     f = [
%         % Position derivatives
%         vx;
%         vy;
%         vz;
% 
%         % Velocity derivatives
%         (T/m) * (c_phi * s_theta * c_psi + s_phi * s_psi);
%         (T/m) * (c_phi * s_theta * s_psi - s_phi * c_psi);
%         (T/m) * c_phi * c_theta - g;
% 
%         % Euler angle derivatives
%         p + q * s_phi * t_theta + r * c_phi * t_theta;
%         q * c_phi - r * s_phi;
%         q * s_phi * sec_theta + r * c_phi * sec_theta;
% 
%         % Angular velocity derivatives
%         (Mp + (Jyy - Jzz) * q * r) / Jxx;
%         (Mq + (Jzz - Jxx) * p * r) / Jyy;
%         (Mr + (Jxx - Jyy) * p * q) / Jzz
%     ];
% end


% function xdot = quadrotor_dynamics_quat(x, u, m, g, J)
%     % State variables
%     p = x(1:3); % Position
%     v = x(4:6); % Velocity
%     q = x(7:10); % Quaternion
%     omega = x(11:13); % Angular velocity
% 
%     % Control inputs
%     T = u(1); % Thrust [N]
%     tau = u(2:4); % Torques [N-m]
% 
%     % Quaternion components
%     q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
% 
%     % Define the symbolic rotation matrix R(q)
%     R = [1 - 2*(q3^2 + q2^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
%          2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
%          2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
% 
% 
%     % Translational dynamics
%     v_dot = [0; 0; -g] + (R * [0; 0; T]) / m;
% 
%     % Rotational dynamics
%     omega_dot = J \ (tau - cross(omega, J * omega));
% 
%     % Quaternion derivative
%     q_dot = 0.5 * quatmultiply(q', [0; omega]')';
% 
%     % Return full state derivative
%     xdot = [v; v_dot; q_dot; omega_dot];
% end
