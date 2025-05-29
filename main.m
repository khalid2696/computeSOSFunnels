clc; clearvars; close all;

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final pose: position and Euler angles (roll, pitch, yaw) in radians

initialPose = [0; 0; 2; 0; 0; 0];   % initial state: origin at height of 2m with zero attitude
finalPose   = [0; -3; 2; 0; 0; 0];   % desired final pose

%% Specify Quadrotor Parameters (not defining these will result in an error)
quadParameters.m = 0.5;        % mass (kg)
quadParameters.g = 9.81;       % gravity (m/s^2)
quadParameters.J = [0.01, 0.01, 0.018]; % moment of inertia (kg⋅m^2) %[4.856e-3, 4.856e-3, 8.801e-3]

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 5;
numTimeSteps = 25;         % number of time samples

drawFlag = 1; % 1: if you want to plot results, 0: otherwise
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

%% for debugging
%load('./precomputedData/nominalTrajectory.mat');
%time_instances(end)
%x_nom(:,1)'
%x_nom(:,end)'

keyboard
%% Design a time-varying LQR feedback controller

%Cost matrices
Q = 0.5*diag([10, 10, 1]); % State cost
R = 0.5*diag([1, 0.1]); % Control cost
%optionally specify terminal cost matrix scaling, P_f = terminalRegionScaling*P_LQR (infinite-time LQR at final state)
terminalRegionScaling = 10; % Terminal constraint cost
                            %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            %and hence increase the terminal cost term (improve tracking/convergence performance) 
                            % most probably, values greater than 1 would work

run("Step2_FeedbackControllerSythesis.m");
disp('- - - - - - -'); disp(" ");

% for debugging
load('./precomputedData/LQRGainsAndCostMatrices.mat');
keyboard

%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

clc; close all;
startTimeIndex = 1; %start time for the rollouts
%startMaxPerturbation = 0.2; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
run("./utils/checkClosedLoop_MCRollouts.m");

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

run("Step3_getDeviationDynamics.m");
%Note: comment out lines 118-122 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% [Optional] Load all the saved files for further analysis

clearvars; close all; 
addpath('./lib/');

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');
%load('./precomputedData/deviationDynamics.mat');

plotOneLevelSet_2D(x_nom, P);
plotOneLevelSet_3D(x_nom, P);

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