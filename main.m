clc; clearvars; close all;

%keyboard
%parpool; %initialise parallel processing %if clusters available and toolboox installed

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

%% [INPUT] specify initial and final state: [pos, vel, theta, omega]
% Convention: theta = 0 -- vertically down (stable), theta = pi -- vertically up (unstable) 


initialState = [0; 0; 0; 0;];  % initial state: at origin, vertically down
finalState   = [3; 0; 0; 0;];  % desired final state

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

maxTimeHorizon = 20; %10 for swing down
numTimeSteps = 25;  % number of time samples %N = 25 for top and down balance
                     %N = 100 for swing-up and swing-down

drawFlag = 1; % 1: if you want to plot results, 0: otherwise
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

% for debugging
load('./precomputedData/nominalTrajectory.mat');
x_nom(:,1)'
x_nom(:,end)'

%% for debugging purposes
%run('./utils/checkOpenLoop_IVP.m');
%run('./utils/checkFeedbackController.m');

%keyboard;

%% Design a time-varying LQR feedback controller

close all
%upsamplingFactor = 20; %finer discretization to prevent integration error build-up
                        %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)

%Cost matrices
% State order: [px; vx; theta; omega]

% Q = diag([5, 0.1, 10, 0.1]);
% R = 1000; %0.1
% %P_f = Q;
% terminalRegionScaling = 100;

Q = diag([5, 0.1, 10, 0.1]);
R = 0.1; %0.1 for down-balance (theta = 0), and %10 for top-balance (theta = pi)
terminalRegionScaling = 1;


run("Step2_FeedbackControllerSynthesis.m");
disp('- - - - - - -'); disp(" ");

% for debugging
load('./precomputedData/LQRGainsAndCostMatrices.mat');

% %condition number (kappa) < 1e2 or 1e3 is good (k=1 is perfect condition and k > 1e3 is ill-conditioned)
% for k=1:1:size(P,3)
%     matrix_condition_number(P(:,:,k))
% end
% 
% keyboard;
%% Additionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all; clearvars;

numSamples = 250;
startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 1e-1; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
upsamplingFactor = 1; %finer discretization to prevent integration error build-up
                       %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)

run("./utils/checkClosedLoop_MCRollouts.m");
drawnow

% for k = 1:50
% 
%     startTimeIndex = k
%     startMaxPerturbation = 0.05;
% 
%     clearvars -except startMaxPerturbation startTimeIndex
% 
%     run("./utils/checkClosedLoop_MCRollouts.m");
%     drawnow
% end
% 
keyboard
%% [Optional] Load all the saved files for further analysis

% clearvars; close all;
% addpath('./lib/');
% 
% load('./precomputedData/nominalTrajectory.mat');
% load('./precomputedData/LQRGainsAndCostMatrices.mat');
% 
% %following functions can take in an additional (optional) argument
% %specifiying the projection dimensions - for example: [1 2], [1 3], etc.
% projectionDims_2D = [1 3];
% plotOneLevelSet_2D(x_nom, P, projectionDims_2D); 
% axis normal
% 
% plotOneLevelSet_3D(x_nom, P);
% daspect([1 1 0.5]);

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

clearvars

order = 3; %order of Taylor expansion
run("Step3_getDeviationDynamics.m");
%Note: comment out lines 109-112 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% Time-conditioned invariant set analysis (with temporal-dependance)
% disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
% clearvars; close all; 
% 
% drawFlag = 1;
% 
% % Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% % Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
% %             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
% %             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
% rhoInitialGuessConstant = 1e-5; %[TUNEABLE] rho_0: decrease value if initial guess fails, 
%                                 % keep it greater than 0!
% rhoInitialGuessExpCoeff = -2; %[TUNEABLE] c: increase value if initial guess fails
% 
% %[0.3, 0.5] %for top balance and down balance
% 
% usageMode = 'feasibilityCheck'; %run just for an initial feasibility check
% try                               
%     run("Step4_computeTimeSampledInvarianceCertificates.m");
% catch
%     disp('Could not find a successful initial guess to start the alternation scheme!');
% end
% 
% keyboard
%% optimize once the feasibility check passes through

disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
clearvars; close all; 

drawFlag = 1;

%specify SOS program hyperparameters
maxIter = 1; %maximum number of iterations
% Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
%             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
%             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
rhoInitialGuessConstant = 0.2; % [TUNEABLE] rho_0: decrease value if initial guess fails, 
                               % keep it greater than 0!
rhoInitialGuessExpCoeff = 0; % [TUNEABLE] c: increase value if initial guess fails

%rho0 = 0.2 and c = 0 seems to give good initial feasible guess (top balance, theta = pi)
%for Q = diag([5, 0.1, 10, 0.1]); R = 10; terminalRegionScaling = 1;

%rho0 = 0.3 and c = 0.5 seems to give good initial feasible guess (down balance, theta = 0)
%for Q = diag([5, 0.1, 10, 0.1]); R = 0.1; terminalRegionScaling = 1;


usageMode = 'shapeOptimisation'; %will have to workshop a name for this!
run("Step4_computeTimeSampledInvarianceCertificates.m");

disp('- - - - - - -'); disp(" ");

%% [Optional] Plot computed funnels
close all

projectionDims_2D = [1 3]; projectionDims_3D = [1 2 4];
run("./utils/plottingScript.m");

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
    omega = x(4);
    
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
        (F + m*L*omega^2*s_theta + m*g*s_theta*c_theta) / denom;
        omega;
        (-F*c_theta - m*L*omega^2*s_theta*c_theta - (M + m)*g*s_theta) / (L * denom)
    ];
end

% Plotting functions
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