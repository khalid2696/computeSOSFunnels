%clc; clearvars; close all

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

% Cost matrices

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% two options for terminal cost:
% 1 - Fixed, input matrix --> define Pf in main file
% 2 - From time-invariant LQR --> computed downstream in this script

% Pf = Q; %Option 1 --> not preferred (requires good sense of the system
%dynamics and scalings between the various constituents of the state vector)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

if ~exist('Q','var')
    Q = eye(size(x_nom,1)); % State cost
end

if ~exist('R','var')
    R = eye(size(u_nom,1)); % Control cost
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

f = dynamicsFnHandle(x, u);

if length(f) ~= length(x)
    error('State vector (x) and dynamics fn (f) are of not same length!');
end

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

%clearvars;
%% Function definitions

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