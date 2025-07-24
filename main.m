clc; clearvars; close all;

%parpool; %initialise parallel processing %if clusters available and toolboox installed

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final state: [pos, vel, theta, omega]
% Convention: theta = 0 -- vertically down (stable), theta = pi -- vertically up (unstable) 

initialState = [0; 0; pi; 0;];   % initial state: at origin, vertically down
finalState   = [2; 0; pi; 0;];   % desired final state

%modify lines 65-67 of ./lib/getNominalTrajectory_using_DirectCollocation.m to impose 
%theta constraints accordingly (based on whether it's upright or hanging down)

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

%% Design a time-varying LQR feedback controller

%upsamplingFactor = 20; %finer discretization to prevent integration error build-up
                        %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)

%Cost matrices
% State order: [px; vx; theta; omega]
Q = diag([5, 0.1, 10, 0.1]);
R = 0.1;
terminalRegionScaling = 1;

run("Step2_FeedbackControllerSynthesis.m");
disp('- - - - - - -'); disp(" ");

% for debugging
%load('./precomputedData/LQRGainsAndCostMatrices.mat');

keyboard;
%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all;
startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 0.1; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
run("./utils/checkClosedLoop_MCRollouts.m");

keyboard
%% [Optional] Load all the saved files for further analysis

clearvars; close all;
addpath('./lib/');

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');

%following functions can take in an additional (optional) argument
%specifiying the projection dimensions - for example: [1 2], [1 3], etc.
projectionDims_2D = [1 3];
plotOneLevelSet_2D(x_nom, P, projectionDims_2D); 
axis normal

plotOneLevelSet_3D(x_nom, P);
daspect([1 1 0.5]);

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

%order = 3; %order of Taylor expansion
run("Step3_getDeviationDynamics.m");
%Note: comment out lines 109-112 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% Time-conditioned invariant set analysis (with temporal-dependance)
disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
clearvars; close all; 

%specify SOS program hyperparameters
maxIter = 1; %maximum number of iterations

% Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
%             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
%             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
rhoInitialGuessConstant = 0.01; %[TUNEABLE] rho_0: decrease value if initial guess fails, 
                                % keep it greater than 0!
rhoInitialGuessExpCoeff = -0.05; %[TUNEABLE] c: increase value if initial guess fails

usageMode = 'feasibilityCheck'; %run just for an initial feasibility check
try                               
    run("Step4_computeTimeSampledInvarianceCertificates.m");
catch
    disp('Could not find a successful initial guess to start the alternation scheme!');
end

%optimize once the feasibility check passes through
usageMode = 'shapeOptimisation'; %will have to workshop a name for this!
run("Step4_computeTimeSampledInvarianceCertificates.m");

disp('- - - - - - -'); disp(" ");

%% [Optional] Plot computed funnels
close all

projectionDims_2D = [1 3]; projectionDims_3D = [1 2 4];
run("./utils/plottingScript.m");

%% [Optional] Verify the theoretical bounds (from SOS programming) with empirical bounds (using MC rollouts) 

startTimeIndex = 1;
numSamples = 100;

run("./utils/empiricallyVerifyInvariance.m");

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
    %denom = M + m*(1 - c_theta^2);
    denom = M + m*s_theta^2;
    
    % State derivatives
    f = [
        v_x;
        (F + m*L*theta_dot^2*s_theta - m*g*s_theta*c_theta) / denom;
        theta_dot;
        (-F*c_theta - m*L*theta_dot^2*s_theta*c_theta + (M + m)*g*s_theta) / (L * denom)
    ];
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