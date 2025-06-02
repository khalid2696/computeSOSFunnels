clc; clearvars; close all;

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final state: [pos, vel, theta, omega]
% Convention: theta = 0 -- vertically down (stable), theta = pi -- vertically up (unstable) 

initialState = [0; 0; 0; 0;];   % initial state: at origin, vertically down
finalState   = [1; 0; pi; 0;];   % desired final state

%% Specify Cart-Pole Parameters (not defining these will result in an error)

cartPoleParameters.M = 1.0;    % cart mass (kg)
cartPoleParameters.m = 0.1;    % pendulum mass (kg)
cartPoleParameters.L = 0.5;    % pendulum length (m)
cartPoleParameters.g = 9.81;   % gravity (m/s^2)

%% Define the system dynamics as a function handle
dynamicsFnHandle = @(x, u) cartpole_dynamics(x, u, cartPoleParameters);

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 10;
numTimeSteps = 50;         % number of time samples

drawFlag = 1; % 1: if you want to plot results, 0: otherwise
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

% for debugging
load('./precomputedData/nominalTrajectory.mat');
x_nom(:,1)'
x_nom(:,end)'

keyboard;
%% Design a time-varying LQR feedback controller

%Cost matrices
% State order: [px; py; pz; vx; vy; vz; phi; theta; psi; p; q; r]
%Q = ;
%R = ;
%terminalRegionScaling = 10;
run("Step2_FeedbackControllerSythesis.m");
disp('- - - - - - -'); disp(" ");

% for debugging
%load('./precomputedData/LQRGainsAndCostMatrices.mat');

%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all;
startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 0.1; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
run("./utils/checkClosedLoop_MCRollouts.m");

%% [Optional] Load all the saved files for further analysis

clearvars; close all;
addpath('./lib/');

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');

plotOneLevelSet_2D(x_nom, P);
%plotOneLevelSet_3D(x_nom, P);

keyboard;
%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

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

% Define the system dynamics: cartpole
% x = [p_x; v_x; theta; theta_dot]
% u = F
function f = cartpole_dynamics(x, u, cartPoleParameters)
    % Numerical evaluation of cartpole dynamics
    
    
    %extract parameters
    M = cartPoleParameters.M; m = cartPoleParameters.m;
    L = cartPoleParameters.L; g = cartPoleParameters.g;

    % Extract states
    p_x = x(1);
    v_x = x(2);
    theta = x(3);
    theta_dot = x(4);
    
    % Control input
    F = u;
    
    % Define trigonometric functions
    s_theta = sin(theta);
    c_theta = cos(theta);
    
    % Common denominator
    denom = M + m - m * c_theta^2;
    
    % State derivatives
    f = [
        v_x;
        (F + m * L * theta_dot^2 * s_theta - m * g * s_theta * c_theta) / denom;
        theta_dot;
        (-F * c_theta - m * L * theta_dot^2 * s_theta * c_theta + (M + m) * g * s_theta) / (L * denom)
    ];
end

function plotOneLevelSet_2D(x_nom, ellipsoidMatrix)
    figure; hold on; grid on; %axis equal;

    P = plottingFnsClass();
    
    projectionDims = [1 3];

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
    ylabel('\theta');
end

function plotOneLevelSet_3D(x_nom, ellipsoidMatrix)
    figure; view(3);
    hold on; grid on; %axis equal;

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
    ylabel('v_x');
    zlabel('\theta');
end