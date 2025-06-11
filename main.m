clc; clearvars; close all;
    
%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final pose: position and Euler angles (roll, pitch, yaw) in radians

initialPose = [0; 0; 2; 0; 0; 0];   % initial state: origin at height of 2m with zero attitude
finalPose   = [2; 4; 2; 0; 0; 0];   % desired final pose

%% Specify Quadrotor Parameters (not defining these will result in an error)
quadParameters.m = 0.5;        % mass (kg)
quadParameters.g = 9.81;       % gravity (m/s^2)
quadParameters.J = [0.01, 0.01, 0.018]; % moment of inertia (kg⋅m^2) %[4.856e-3, 4.856e-3, 8.801e-3]

%% Define the system dynamics as a function handle
dynamicsFnHandle = @(x, u) quadrotor_dynamics(x, u, quadParameters);

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 5;
numTimeSteps = 25;         % number of time samples

drawFlag = 1; % 1: if you want to plot results, 0: otherwise
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

%% Design a time-varying LQR feedback controller

%Cost matrices
% State order: [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
Q = 0.1*diag([
    50,  50,  50, ...   % Position weights [px, py, pz] - high for tracking
    5,   5,   5,  ...   % Velocity weights [vx, vy, vz] - moderate
    10,   10,   10,  ...   % Attitude weights [phi, theta, psi] - roll/pitch higher than yaw
    1,    1,    1  ]);     % Angular rate weights [p, q, r] - low, let inner loop handle

% Q_aggressive = diag([200, 200, 200, 20, 20, 20, 80, 80, 10, 2, 2, 2]);
% R_conservative = diag([0.05, 20, 20, 30]);

% Control order: [T; Mx; My; Mz]
R = diag([
    1, ...    % Thrust - relatively low, thrust is "cheap"
    30,  ...    % Roll moment Mx - higher, moments are more "expensive"
    30,  ...    % Pitch moment My
    50  ]);     % Yaw moment Mz - highest, yaw control typically less aggressive

%optionally specify terminal cost matrix scaling, P_f = terminalRegionScaling*P_LQR (infinite-time LQR at final state)
terminalRegionScaling = 20; % Terminal constraint cost
                            %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            %and hence increase the terminal cost term (improve tracking/convergence performance) 
                            % most probably, values greater than 1 would work

% Tuning guidelines:
% For better position tracking:
% 
% Increase position weights (px, py, pz): 100 to 200
% Increase velocity weights if oscillations occur: 10 to 20
% 
% For smoother control:
% 
% Increase R values (especially moments): Mx, My from 10 to 50
% Decrease attitude weights if too aggressive: 50 to 20
% 
% Physical reasoning:
% 
% Position weights (100): High because position tracking is primary objective
% Attitude weights (50, 50, 20): Roll/pitch couple to xy-position, so weighted higher than yaw
% Thrust weight (0.1): Low because thrust changes are energetically cheap
% Moment weights (10, 10, 20): Higher because moments require more energy and we want smooth attitude control

run("Step2_FeedbackControllerSynthesis.m");
disp('- - - - - - -'); disp(" ");

% for debugging
load('./precomputedData/LQRGainsAndCostMatrices.mat');

keyboard;
%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all;
startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 1; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
run("./utils/checkClosedLoop_MCRollouts.m");

%% [Optional] Load all the saved files for further analysis

clearvars; %close all;
addpath('./lib/');

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');

plotOneLevelSet_2D(x_nom, P);
plotOneLevelSet_3D(x_nom, P);

keyboard;
%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

order = 3; %order of Taylor expansion
run("Step3_getDeviationDynamics.m");
%Note: comment out lines 109-112 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% Time-conditioned invariant set analysis (with time-dependance)
disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
clearvars; close all; 
keyboard;

%specify SOS program hyperparameters
maxIter = 1; %maximum number of iterations
rhoGuessChoice = 'const'; %options: 'const' [DEFAULT] and 'exp'
                          %[USE 'exp' ONLY IF 'const' IS INFEASIBLE]
% Option1: constant rho_guess
rhoInitialGuessConstant = 0.4; %[TUNEABLE] decrease value if initial guess fails, 
                               %keep value between 0 and 1

% Option2: Exponentially (quickly) increasing rho_guess 
rhoInitialGuessExpCoeff = 1.5; %[TUNEABLE] increase value if initial guess fails
                               %keep value greater than 0 (increasing fn)

run("Step4_computeTimeSampledInvarianceCertificates.m");
disp('- - - - - - -'); disp(" ");

%% [Optional] Plot computed funnels
run("./utils/plottingScript.m");

%% [Optional] Verify the theoretical bounds (from SOS programming) with empirical bounds (using MC rollouts) 

%startTimeIndex = 1;
%numSamples = 100;

%run("./utils/empiricallyVerifyInvariance.m");

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
    t_theta = tan(theta);
    sec_theta = sec(theta);
    
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

function plotOneLevelSet_2D(x_nom, ellipsoidMatrix)
    figure; hold on; grid on; axis equal;

    P = plottingFnsClass();
    
    projectionDims = [1 2];

    %plot ellipsoidal invariant sets in 2D
    for k=1:1:length(x_nom)
        M = ellipsoidMatrix(:,:,k);
        M_xy = P.project_ellipsoid_matrix_2D(M, projectionDims);
        %center = x_nom(:,k);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k)]';
        %plotEllipse(center, M_xy);
        P.plotEllipse(center, M_xy);
    end 
    
    %nominal trajectory
    plot(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),'--b');

    %formatting
    title('1-level set of cost-to-go matrices P, along the nominal trajectory');
    xlabel('p_x');
    ylabel('p_y');
end

function plotOneLevelSet_3D(x_nom, ellipsoidMatrix)
    figure; view(3);
    hold on; grid on; axis equal;

    P = plottingFnsClass();

    projectionDims = [1 2 3];
    
    %plot ellipsoidal invariant sets in 3D
    for k=1:1:length(x_nom)
        M = ellipsoidMatrix(:,:,k);
        M_xyz = P.project_ellipsoid_matrix_3D(M, projectionDims);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k), x_nom(projectionDims(3),k)]';
        P.plotEllipsoid(center, M_xyz);
    end 
    
    %nominal trajectory
    plot3(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),x_nom(projectionDims(3),:),'--b');

    %formatting
    title('1-level set of cost-to-go matrices P, along the nominal trajectory');
    xlabel('p_x');
    ylabel('p_y');
    zlabel('p_z');
end