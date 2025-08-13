clc; clearvars; close all;

%keyboard
%parpool; %initialise parallel processing %if clusters available and toolboox installed

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final pose: position and Euler angles (roll, pitch, yaw) in radians

initialPose = [0; 0; 2; 0; 0; 0];   % initial state: origin at height of 2m with zero attitude
finalPose   = [-2; -5; 2; 0; 0; 0];   % desired final pose

%% Specify Quadrotor Parameters (not defining these will result in an error)
quadParameters.m = 0.7;                  % mass (kg)
quadParameters.g = 9.81;                 % gravity (m/s^2)
quadParameters.J = [2e-3, 2e-3, 3.5e-3]; % moment of inertia (kg⋅m^2) %[4.856e-3, 4.856e-3, 8.801e-3]

%% Define the system dynamics as a function handle
dynamicsFnHandle = @(x, u) quadrotor_dynamics(x, u, quadParameters);

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 15;
numTimeSteps   = 25;         % number of time samples

drawFlag = 1; % 1: if you want to plot results, 0: otherwise
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

%keyboard

%% Design a time-varying LQR feedback controller
clearvars
%upsamplingFactor = 1; %finer discretization to prevent integration error build-up
                       %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)

%Cost matrices
% State order: [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
% Q = 0.01*diag([
%     50,  50,  80, ...   % Position weights [px, py, pz]
%     5,   5,   8,  ...   % Velocity weights [vx, vy, vz]
%     20,  20,  10,  ...  % Attitude weights [phi, theta, psi] - roll/pitch higher than yaw
%     100, 100,  100  ]);   % Angular rate weights [p, q, r]

% Q = 0.1*diag([
%     100,  100,  400, ...   % Position weights [px, py, pz]
%     4,  4,  25,  ...   % Velocity weights [vx, vy, vz]
%     200,  200,  100,  ...  % Attitude weights [phi, theta, psi] - roll/pitch higher than yaw
%     4, 4,  4  ]);   % Angular rate weights [p, q, r]

Q = 0.1*eye(12);
Q(3,3) = 0.4;
Q(7,7) = 10; Q(8,8) = 10; Q(9,9) = 10;
Q(10,10) = 1; Q(11,11) = 1; Q(12,12) = 1;

% R = 10*diag([
%     0.2, ...    % Thrust
%     1,  ...   % Roll moment Mx
%     1,  ...   % Pitch moment My
%     1  ]);   % Yaw moment Mz

R = 10*eye(4);
R(1,1) = 2;

%optionally specify terminal cost matrix scaling, P_f = terminalRegionScaling*P_LQR (infinite-time LQR at final state)
terminalRegionScaling = 1; % Terminal constraint cost
                            %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            %and hence increase the terminal cost term (improve tracking/convergence performance) 
                            % most probably, values greater than 1 would work

run("Step2_FeedbackControllerSynthesis.m");
disp('- - - - - - -'); disp(" ");

load('./precomputedData/LQRGainsAndCostMatrices.mat');

%condition number (kappa) < 1e2 or 1e3 is good (k=1 is perfect condition and k > 1e3 is ill-conditioned)
for k=1:1:size(P,3)
    matrix_condition_number(P(:,:,k))
end

%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all; clearvars;

startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 0.05; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
upsamplingFactor = 40; %finer discretization to prevent integration error build-up
                       %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)
                         
run("./utils/checkClosedLoop_MCRollouts.m");
drawnow

% for k = 1:25
% 
%     startTimeIndex = k
%     startMaxPerturbation = 0.05;
% 
%     clearvars -except startMaxPerturbation startTimeIndex
% 
%     run("./utils/checkClosedLoop_MCRollouts.m");
%     drawnow
% end

%keyboard
%% [Optional] Load all the saved files for further analysis

% clearvars; %close all;
% addpath('./lib/');
% 
% load('./precomputedData/nominalTrajectory.mat');
% load('./precomputedData/LQRGainsAndCostMatrices.mat');
% 
% %following functions can take in an additional (optional) argument
% %specifiying the projection dimensions - for example: [1 2], [1 3], etc.
% projectionDims_2D = [1 2];
% plotOneLevelSet_2D(x_nom, P, projectionDims_2D); %axis auto or %axis normal
% 
% projectionDims_3D = [1 2 3];
% plotOneLevelSet_3D(x_nom, P, projectionDims_3D);
% 
% projectionDims_3D = [4 5 6];
% plotOneLevelSet_3D(x_nom, P, projectionDims_3D);
% 
% projectionDims_3D = [7 8 9];
% plotOneLevelSet_3D(x_nom, P, projectionDims_3D);
% 
% projectionDims_3D = [10 11 12];
% plotOneLevelSet_3D(x_nom, P, projectionDims_3D);

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

order = 3; %order of Taylor expansion
run("Step3_getDeviationDynamics.m");
%Note: comment out lines 109-112 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% Time-conditioned invariant set analysis (with temporal-dependance)
disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
clc; clearvars; close all; 

drawFlag = 1;

%specify SOS program hyperparameters
maxIter = 1; %maximum number of iterations

% Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
%             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
%             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
rhoInitialGuessConstant = 1e-3; %[TUNEABLE] rho_0: decrease value if initial guess fails, 
                                % keep it greater than 0!
rhoInitialGuessExpCoeff = 3; %[TUNEABLE] c: decrease value if initial guess fails

%1e-3 and 3

usageMode = 'feasibilityCheck'; %run just for an initial feasibility check
try                               
    run("Step4_computeTimeSampledInvarianceCertificates.m");
catch
    disp('Could not find a successful initial guess to start the alternation scheme!');
end

%return

%% optimize once the feasibility check passes through
close all
%drawFlag = 1;
tic
usageMode = 'shapeOptimisation'; %will have to workshop a name for this!
run("Step4_computeTimeSampledInvarianceCertificates.m");
toc

%% [Optional] Plot computed funnels
clc; clearvars; close all

projectionDims_2D = [1 2]; projectionDims_3D = [1 2 3];
run("./utils/plottingScript.m");

%% [Optional] Verify the theoretical bounds (from SOS programming) with empirical bounds (using MC rollouts) 

% startTimeIndex = 1;
% numSamples = 100;
% 
% run("./utils/empiricallyVerifyInvariance.m");

%% Small debug snippets

% alpha = 1e-3;  %default: 0.001
% expCoeff = 1; %default: 0.1
% 
% N = length(time_instances);
% t0 = time_instances(1); tf = time_instances(N);
% rhoScaleIncrements = NaN(size(time_instances));
% for k = 1:N
%     tk = time_instances(k);
%     rhoScaleIncrements(k) = 1 + alpha * exp(expCoeff*(tk - tf)/(t0 - tf));
% end
% 
% plot(time_instances, rhoScaleIncrements); xlabel('time'); ylabel('\rho scaling')
% drawnow;

%% Function definitions

% Define the system dynamics: quadrotor model
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
    % Input: u = [T; Mx; My; Mz]
    T = u(1); Mp = u(2); Mq = u(3); Mr = u(4);
    
    % Trigonometric shortcuts
    c_phi = cos(phi); s_phi = sin(phi);
    c_theta = cos(theta); s_theta = sin(theta);
    c_psi = cos(psi); s_psi = sin(psi);
    t_theta = s_theta/c_theta;
    sec_theta = 1/c_theta;
    
    % Quadrotor dynamics 
    % Assumptions: no aerodynamic drag and gyroscopic coupling due to rotor inertia)
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

% Computes the condition number of a matrix to check for numerical ill-conditioning:
% kappa < 1e2 or 1e3 is good (k=1 is perfect condition and k > 1e3 is ill-conditioned)
function kappa = matrix_condition_number(A, norm_type)
% Input:
%   A - Input matrix (must be square and non-singular)
%   norm_type - (optional) Type of norm to use:
%               '1' or 1 for 1-norm
%               '2' or 2 for 2-norm (default)
%               'inf' for infinity-norm
%               'fro' for Frobenius norm
%
% Output:
%   kappa - Condition number of matrix A

    % Set default norm type if not provided
    if nargin < 2
        norm_type = 2;
    end
    
    % Check if matrix is square
    [m, n] = size(A);
    if m ~= n
        disp('Matrix must be square');
        return
    end
    
    % Check if matrix is singular
    if abs(det(A)) < eps
        warning('Matrix is close to singular. Condition number may be very large.');
    end
    
    % Compute condition number using built-in function
    kappa = cond(A, norm_type);
    
    % Alternative manual calculation (commented out):
    % kappa = norm(A, norm_type) * norm(inv(A), norm_type);
    
end

function plotOneLevelSet_2D(x_nom, ellipsoidMatrix, projectionDims)
    figure; hold on; grid on; axis equal;

    P = plottingFnsClass();
    
    if nargin < 3
        projectionDims = [1 2]; %if not specified, by default x-y projection
    end

    %plot ellipsoidal invariant sets in 2D
    for k=1:1:size(x_nom,2)
        M = ellipsoidMatrix(:,:,k);
        M_xy = P.project_ellipsoid_matrix_2D(M, projectionDims);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k)]';
        P.plotEllipse(center, M_xy);
    end 
    
    %nominal trajectory
    plot(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),'--b');

    %formatting
    title('1-level set of cost-to-go matrices P, along the nominal trajectory');
    xlabel(['x_{', num2str(projectionDims(1)), '}'])
    ylabel(['x_{', num2str(projectionDims(2)), '}'])
end

function plotOneLevelSet_3D(x_nom, ellipsoidMatrix, projectionDims)
    figure; view(3);
    hold on; grid on; axis equal;

    P = plottingFnsClass();

    if nargin < 3
        projectionDims = [1 2 3]; %if not specified, by default x-y-z projection
    end
    
    %plot ellipsoidal invariant sets in 3D
    for k=1:1:size(x_nom,2)
        M = ellipsoidMatrix(:,:,k);
        M_xyz = P.project_ellipsoid_matrix_3D(M, projectionDims);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k), x_nom(projectionDims(3),k)]';
        P.plotEllipsoid(center, M_xyz);
    end 
    
    %nominal trajectory
    plot3(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),x_nom(projectionDims(3),:),'--b');

    %formatting
    title('1-level set of cost-to-go matrices P, along the nominal trajectory');
    xlabel(['x_{', num2str(projectionDims(1)), '}'])
    ylabel(['x_{', num2str(projectionDims(2)), '}'])
    zlabel(['x_{', num2str(projectionDims(3)), '}'])
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