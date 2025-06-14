%clc; clearvars; close all

disp('Synthesizing an LQR controller for tracking the nominal trajectory..');
disp(' ');

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

%finer discretization to prevent integration error build-up
if ~exist('upsamplingFactor','var')
    upsamplingFactor = 10;
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

nx = size(x_nom, 1); nu = size(u_nom, 1); N = length(time_instances);
x = createSymbolicVector('x', nx); %state vector
u = createSymbolicVector('u', nu); %input vector
% [x1, x2, x3] -- [px; py; theta]
% [u1, u2] -- [v, omega]

f = dynamicsFnHandle(x, u);

if length(f) ~= length(x)
    error('State vector (x) and dynamics fn (f) are of not same length!');
end

%% Define terminal cost matrix if not specified using TI-LQR
if ~exist('P_f','var')
    [~, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom(:,end), u_nom(:,end), Q, R);
    P_f = S_f*terminalRegionScaling;
end

%% Interpolate state and input trajectory for a finer discretisation, so that integration errors don't build up

% Number of interpolated points
Nd = N*upsamplingFactor;

% Fine time vector for interpolation - upsampled spacing
t_fine = linspace(time_instances(1), time_instances(end), Nd);

[x_fine, u_fine] = upsample_state_control_trajectories(time_instances, x_nom, u_nom, t_fine);

%% Compute time-sampled TVLQR Gains on a finer discretization

Ts = (t_fine(end)-t_fine(1))/(length(t_fine)-1); 
[K_fine, P_fine] = compute_tvlqr_gains(f, x, u, x_fine, u_fine, Q, R, P_f, Ts);

%% Downsample the trajectories, inputs, cost and gain matrices to match the original time samples

[x_nom, u_nom] = downsample_state_control_trajectories(t_fine, x_fine, u_fine, time_instances);

K = downsample_matrix(K_fine, t_fine, time_instances);
P = downsample_matrix(P_fine, t_fine, time_instances);

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

%Interpolates state and control vectors
function [x_fine, u_fine, t_fine] = upsample_state_control_trajectories(t, x, u, t_fine)

    % Dimensions
    Nd = length(t_fine); n = size(x, 1); m = size(u, 1);

    % Preallocate outputs
    x_fine = zeros(n, Nd);
    u_fine = zeros(m, Nd);

    % Interpolate each row (dimension) of x with cubic spline
    for i = 1:n
        x_fine(i, :) = interp1(t, x(i, :), t_fine, 'spline'); % Cubic interpolation
    end

    % Interpolate each row (dimension) of u with linear interpolation
    for i = 1:m
        u_fine(i, :) = interp1(t, u(i, :), t_fine, 'linear'); % Linear interpolation
    end
end

%downsamples state and control trajectories to match with given (coarse) time samples
function [x_coarse, u_coarse] = downsample_state_control_trajectories(t_fine, x_fine, u_fine, t_coarse)
    
    % Dimensions
    N = length(t_coarse); n = size(x_fine, 1); m = size(u_fine, 1);

    % Downsample to match original time vector t
    x_coarse = zeros(n, N);  % n x N
    u_coarse = zeros(m, N);  % m x N
    
    for i = 1:n
        x_coarse(i, :) = interp1(t_fine, x_fine(i, :), t_coarse, 'spline');  % or 'pchip' if smoother
    end
    
    for i = 1:m
        u_coarse(i, :) = interp1(t_fine, u_fine(i, :), t_coarse, 'linear');
    end
end

%downsamples an input matrix to match with the given (coarse) time samples
function M_coarse = downsample_matrix(M_fine, t_fine, t_coarse)
    % M_fine  : d1 × d2 × Nd
    % t_fine  : 1 × Nd
    % t_coarse: 1 × N (time vector corresponding to downsampling)
    
    [dim1, dim2, ~] = size(M_fine); N = length(t_coarse);
    
    M_coarse = zeros(dim1, dim2, N);
    for i = 1:dim1
        for j = 1:dim2
            M_coarse(i, j, :) = interp1(t_fine, squeeze(M_fine(i, j, :)), t_coarse, 'linear');
            %note: 'linear' and 'spline' (cubic) don't seem to have much difference here
        end
    end
end

%Interpolates an input matrix
function M_fine = upsample_matrix(M_coarse, t_coarse, t_fine)
    % M_coarse: d1 × d2 × N
    % t_coarse: 1 × N
    % t_fine  : 1 × Nd (time vector corresponding to upsampling)
    
    [dim1, dim2, ~] = size(M_coarse); N_fine = length(t_fine);
    
    M_fine = zeros(dim1, dim2, N_fine);
    for i = 1:dim1
        for j = 1:dim2
            M_fine(i, j, :) = interp1(t_coarse, squeeze(M_coarse(i, j, :)), t_fine, 'linear');
        end
    end
end

%% Usage Note:
%Vq = interp1(X,V,Xq) %interpolates to find Vq, the values of the
%                     %underlying function V=F(X) at the query points Xq.