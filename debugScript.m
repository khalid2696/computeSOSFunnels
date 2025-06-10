
% Parameters
cartPoleParameters.M = 1.0;    % cart mass (kg)
cartPoleParameters.m = 0.1;    % pendulum mass (kg)
cartPoleParameters.L = 0.5;    % pendulum length (m)
cartPoleParameters.g = 9.81;   % gravity (m/s^2)

% Returns a numerical function handle for use in ODE solvers
f_handle = @(x, u) cartpole_dynamics_fn(x, u, cartPoleParameters);

% Example: Linearize around upright equilibrium (theta = pi)
x_eq = [0; 0; pi; 0];  % cart at origin, pendulum upright, at rest
u_eq = 0;              % no force needed at equilibrium

% Test the dynamics
xdot_eq = f_handle(x_eq, u_eq);
fprintf('State derivative at upright equilibrium: [%.3f, %.3f, %.3f, %.3f]\n', xdot_eq);

syms x1 x2 x3 x4 u
x = [x1 x2 x3 x4];

f_sym = f_handle(x, u);

% For linearization, you would compute Jacobians:
A = jacobian(f_sym, x) %evaluated at (x_eq, u_eq)
B = jacobian(f_sym, u) %evaluated at (x_eq, u_eq)

fprintf('\nFor TVLQR, compute Jacobians A and B numerically or symbolically\n');
fprintf('around your reference trajectory.\n');

function f = cartpole_dynamics_analytical()
    % Cartpole Dynamics in Analytical Closed Form
    % State: x = [x; x_dot; theta; theta_dot]
    % Input: u = F (horizontal force on cart)
    % Convention: theta = 0 (down), theta = pi (up), positive counterclockwise
    
    % Physical parameters
    M = 1.0;    % cart mass (kg)
    m = 0.1;    % pendulum mass (kg)
    L = 0.5;    % pendulum length to center of mass (m)
    g = 9.81;   % gravity (m/s^2)
    
    % State variables (symbolic for analytical form)
    syms x x_dot theta theta_dot real
    
    % Control input
    syms F real
    
    % Trigonometric functions
    s_theta = sin(theta);
    c_theta = cos(theta);
    
    % Common terms in dynamics
    denominator = M + m - m * c_theta^2;
    
    % Analytical closed-form dynamics f(x,u) = x_dot
    f = [
        % Cart position derivative
        x_dot;
        
        % Cart acceleration (from coupled equations of motion)
        (F + m * L * theta_dot^2 * s_theta - m * g * s_theta * c_theta) / denominator;
        
        % Pendulum angle derivative  
        theta_dot;
        
        % Pendulum angular acceleration (from coupled equations of motion)
        (-F * c_theta - m * L * theta_dot^2 * s_theta * c_theta + (M + m) * g * s_theta) / (L * denominator)
    ];
    
    % Display the analytical expressions
    fprintf('Cartpole Dynamics f(x,u):\n');
    fprintf('x_dot_1 = x_dot\n');
    fprintf('x_dot_2 = (F + m*L*theta_dot^2*sin(theta) - m*g*sin(theta)*cos(theta)) / (M + m - m*cos(theta)^2)\n');
    fprintf('x_dot_3 = theta_dot\n');
    fprintf('x_dot_4 = (-F*cos(theta) - m*L*theta_dot^2*sin(theta)*cos(theta) + (M+m)*g*sin(theta)) / (L*(M + m - m*cos(theta)^2))\n');
    fprintf('\nParameters: M=%.1f kg, m=%.1f kg, L=%.1f m, g=%.2f m/s^2\n', M, m, L, g);
end


