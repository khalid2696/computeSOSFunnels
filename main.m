clc; clearvars; close all;

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% [INPUT] specify initial and final states

xinitial = [0; 0; pi/2];   % initial state: origin, pointing North 
xfinal   = [2; 5; pi/2];   % desired final state 

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 10;
numTimeSteps = 25;         % number of time samples

drawFlag = 1; %uncomment this if you want to plot results
run("Step1_NominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

%% Design a time-varying LQR feedback controller

%Cost matrices
Q = 0.5*diag([10, 10, 1]); % State cost
R = 0.5*diag([1, 0.1]); % Control cost
%optionally specify terminal cost matrix scaling, P_f = terminalRegionScaling*P_LQR (infinite-time LQR at final state)
terminalRegionScaling = 10; % Terminal constraint cost
                            %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            %and hence increase the terminal cost term (improve tracking/convergence performance) 
                            % most probably, values greater than 1 would work

run("Step2_LQR.m");
disp('- - - - - - -'); disp(" ");

%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing
%startTimeIndex = 1; %start time for the rollouts
%initial_state_covariance = (1/3)*(M(:,:,1)/rhoGuess(1))^(-1/2); %specify the covariance for sampling initial states
%run("./utils/checkClosedLoop_usingMCRollouts.m");

% Or, alternatively visualise the flow field of the closed loop system 
%run("./utils/visualise_flowField.m");

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

run("Step3_getDeviationDynamics.m");
%Note: comment out lines 118-122 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% [Optional] Load all the saved files for further analysis

%clearvars;

%load('./precomputedData/nominalTrajectory.mat');
%load('./precomputedData/LQRGainsAndCostMatrices.mat');
%load('./precomputedData/deviationDynamics.mat');

%plotOneLevelSet(x_nom, P);

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

%run("Step4_computeTimeSampledInvarianceCertificates.m");
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

function plotOneLevelSet(x_nom, ellipsoidMatrix)
    figure
    hold on;
    grid on; 
    axis equal;
    
    for k=1:1:length(x_nom)
        M = ellipsoidMatrix(:,:,k);
        M_xy = project_ellipsoid_matrix(M, [1 2]);
        center = x_nom(:,k);
        plotEllipse(center, M_xy)
    end 
    
    title('1-level set of cost-to-go matrices P, along the nominal trajectory');
    xlabel('p_x');
    ylabel('p_y');
    plot(x_nom(1,:),x_nom(2,:),'--b');
end

function M_2d = project_ellipsoid_matrix(M, projection_dims)
    % Input:
    % M: nxn matrix defining the n-dimensional ellipsoid x^T M x < 1
    % projection_dims: 2-element vector specifying which dimensions to project onto
    %                  (e.g., [1 2] for xy-plane, [1 3] for xz-plane)
    
    n = size(M, 1); %get the dimensionality of matrix M

    basisMatrix = zeros(n,2);
    basisMatrix(projection_dims(1),1) = 1;
    basisMatrix(projection_dims(2),2) = 1;

    M_2d = inv(basisMatrix' *inv(M) * basisMatrix);
end

function plotEllipse(center, ellipseMatrix)
    
    %plot an ellipse from which initial states are sampled
    ellipseCenter = center(1:2); % 2D center of the ellipsoid
    [eig_vec, eig_val] = eig(ellipseMatrix);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(ellipseCenter(1) + rotated_ellipse(1, :), ...
         ellipseCenter(2) + rotated_ellipse(2, :), ...
         '-k', 'LineWidth', 1.2);  
end
