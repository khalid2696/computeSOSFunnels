clc; clearvars; close all

% Quadrotor Parameters
quadParameters.m = 0.5;        % mass (kg)
quadParameters.g = 9.81;       % gravity (m/s^2)
quadParameters.J = [0.01, 0.01, 0.018]; % moment of inertia (kg⋅m^2) %[4.856e-3, 4.856e-3, 8.801e-3]

% State variables
syms px py pz vx vy vz phi theta psi p q r real
x = [px py pz vx vy vz phi theta psi p q r]';

% Control inputs
syms T Mp Mq Mr real
u = [T Mp Mq Mr];

%fn handle for the dynamics
quad_dynamics = @(x, u) quadrotor_dynamics(x, u, quadParameters);

f = quadrotor_dynamics(x, u, quadParameters)

% Jacobian with respect to state
A = jacobian(f, x);

% Jacobian with respect to input  
B = jacobian(f, u);


% Quadrotor dynamics 
% Assumptions: no aerodynamic drag and gyroscopic coupling due to rotor inertia)
function f = quadrotor_dynamics(x, u, quadParameters)
    % Numerical version with specific parameter values
    
    %extracting quadrotor model parameters
    m = quadParameters.m; g = quadParameters.g;
    Jxx = quadParameters.J(1); Jyy = quadParameters.J(2); Jzz = quadParameters.J(3);
    
    %assigning state variables for ease of usage
    % State: x = [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
    px = x(1); py = x(2); pz = x(3);
    vx = x(4); vy = x(5); vz = x(6);
    phi = x(7); theta = x(8); psi = x(9);
    p = x(10); q = x(11); r = x(12);

    %assigning input variables for ease of usage
    % Input: u = [T; Mx; My; Mz]
    T = u(1); Mp = u(2); Mq = u(3); Mr = u(4);
    
    % Trigonometric shortcuts
    c_phi = cos(phi); s_phi = sin(phi);
    c_theta = cos(theta); s_theta = sin(theta);
    c_psi = cos(psi); s_psi = sin(psi);
    t_theta = tan(theta);
    sec_theta = sec(theta);
    
    % Closed-form dynamics
    f = [
        % Position derivatives
        vx;
        vy;
        vz;
        
        % Velocity derivatives
        (T/m) * (c_phi * s_theta * c_psi + s_phi * s_psi);
        (T/m) * (c_phi * s_theta * s_psi - s_phi * c_psi);
        (T/m) * c_phi * c_theta - g;
        
        % Euler angle derivatives
        p + q * s_phi * t_theta + r * c_phi * t_theta;
        q * c_phi - r * s_phi;
        q * s_phi * sec_theta + r * c_phi * sec_theta;
        
        % Angular velocity derivatives
        (Mp + (Jyy - Jzz) * q * r) / Jxx;
        (Mq + (Jzz - Jxx) * p * r) / Jyy;
        (Mr + (Jxx - Jyy) * p * q) / Jzz
    ];
end