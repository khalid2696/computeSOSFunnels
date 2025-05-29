clc; clearvars; close all

%% Add directories
addpath('./lib/');

%% Nominal trajectory and nominal/open-loop control input (feedforward term)
load('./precomputedData/nominalTrajectory.mat');

% -- Status message: Quantities at our disposal now -- %

% time_instances   : time horizon (sampled)        : 1 x N
% x_nom            : nominal state trajectory      : n_x x N
% u_nom            : feedforward input tape        : n_u x N
% dynamicsFnHandle : the function handle of the system dynamics

%% Specify parameters or Inherit them if they exist in the wrapper file

if ~exist('quadParameters','var')
    quadParameters.m = 0.5;        % mass (kg)
    quadParameters.g = 9.81;       % gravity (m/s^2)
    quadParameters.J = [0.01, 0.01, 0.018]; % moment of inertia (kg⋅m^2) 
                       %[4.856e-3, 4.856e-3, 8.801e-3]
end

% Cost matrices

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% two options for terminal cost:
% 1 - Fixed, input matrix --> define Pf in main file
% 2 - From time-invariant LQR --> computed downstream in this script

% Pf = Q; %Option 1 --> not preferred (requires good sense of the system
%dynamics and scalings between the various constituents of the state vector)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

if ~exist('Q','var')
    Q = eye(12); % State cost
end

if ~exist('R','var')
    R = eye(4); % Control cost
end

if ~exist('terminalRegionScaling','var')
    terminalRegionScaling = 10; % Terminal constraint cost
     %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            % most probably, values greater than 1 would work
end

%% Define the symbolic system dynamics - states, input and the nonlinear function

nx = size(x_nom, 1); nu = size(u_nom, 1);
x = createSymbolicVector('x', nx); %state vector
u = createSymbolicVector('u', nu); %input vector
% [x1, x2, .., x12] -- [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
% [u1, u2, u3, u4]  -- [T; Mp; Mq; Mr]

f = quadrotor_dynamics(x, u, quadParameters);

%% Define terminal cost matrix if not specified using TI-LQR
if ~exist('P_f','var')
    [~, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom(:,end), u_nom(:,end), Q, R);
    P_f = S_f*terminalRegionScaling;
end

%% Compute time-sampled TVLQR Gains

% Sampling-time from the nominal trajectory
Ts = (time_instances(end)-time_instances(1))/(length(time_instances)-1); 

[K, P] = compute_tvlqr_gains(f, x, u, x_nom, u_nom, Q, R, P_f, Ts);

disp('Finished synthesizing a time-varying LQR stabilizing feedback controller');
disp(' ');

%% save the nominal trajectory and LQR gains and cost-to-go matrices
f_sym = f;
save('./precomputedData/LQRGainsAndCostMatrices.mat', 'time_instances', 'K', 'P' , 'f_sym');

disp('Saved the time-sampled LQR gains and cost-to-go matrices to a file!');
disp(' ');

clearvars;
%% Function definitions

% Quadrotor dynamics 
% Assumptions: no aerodynamic drag and gyroscopic coupling due to rotor inertia)
function f = quadrotor_dynamics(x, u, quadParameters)
    % Numerical version with specific parameter values
    
    %extracting quadrotor model parameters
    m = quadParameters.m; g = quadParameters.g;
    Jxx = quadParameters.J(1); Jyy = quadParameters.J(2); Jzz = quadParameters.J(3);
    
    %assigning state variables for ease of usage
    % State: x = [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
    %px = x(1); py = x(2); pz = x(3);
    vx = x(4); vy = x(5); vz = x(6);
    phi = x(7); theta = x(8); psi = x(9);
    p = x(10); q = x(11); r = x(12);

    %assigning input variables for ease of usage
    % Input: u = [T; Mp; Mq; Mr]
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


%inputs 
% varName: string -- name of vector
% n: int -- length of vector
function sym_vector = createSymbolicVector(vecSymbol, n)
    % Generate variable names dynamically
    varNames = arrayfun(@(i) [vecSymbol, num2str(i)], 1:n, 'UniformOutput', false);
    
    % Define these as symbolic variables
    syms(varNames{:});
    
    % Construct a symbolic vector from these variables
    sym_vector = sym(varNames, 'real');

    sym_vector = sym_vector'; %column matrix by convention -- nx1
end


%an infinite horizon LQR to stabilize the closed loop system at the terminal state
function [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom_f, u_nom_f, Q, R)

    A_symb = jacobian(f, x);
    B_symb = jacobian(f, u);
    
    %Get the linearised matrices (Jacobians) at terminal state
    A_f = double(subs(A_symb, [x; u], [x_nom_f; u_nom_f]));
    B_f = double(subs(B_symb, [x; u], [x_nom_f; u_nom_f]));
    
    [K_f,S_f,~] = lqr(A_f,B_f,Q,R);
end