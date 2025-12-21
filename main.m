clc; clearvars; close all;
%parpool %for parallel processing

%Note: Not defining input-parameters in these files WILL NOT lead to errors 
%      However, the executed scripts will assume some default input values
%WARNING: Will have to use the same variable names as below if we want to
%         pass onto the corresponding scripts

%% Add directories
addpath('./lib/');

<<<<<<< HEAD
%% [INPUT] specify initial and final state: [pos_x, pos_y, theta]
% Convention: theta = pi/2 -- North, 0 -- East

initialState = [ 0; 0; pi/2;];   % initial state: at origin facing North
finalState   = [2; 5; pi/2;];   % desired final state
=======
%% [INPUT] specify   initial and final state: [pos_x, pos_y, theta]
% Convention: theta = pi/2 -- North, 0 -- East

initialState = [ 0; 0; pi/2;];   % initial state: at origin facing North
finalState   = [-2; 4; pi/2;];   % desired final state
>>>>>>> unicycle

%% Define the system dynamics as a function handle
dynamicsFnHandle = @(x, u) unicycle_dynamics(x, u);

%% Compute a nominal trajectory and corresponding feedforward inputs

maxTimeHorizon = 10;
numTimeSteps = 40;         % number of time samples

drawFlag = 0; %uncomment this if you want to plot results
run("Step1_computeNominalTrajectory.m");
disp('- - - - - - -'); disp(" ");

%keyboard
%% Design a time-varying LQR feedback controller

upsamplingFactor = 10; %finer discretization to prevent integration error build-up
                       %finer num of samples = upsamplingFactor*numTimeSamples (temporarily)

%Cost matrices
Q = 0.1*diag([10, 10, 1]); % State cost
R = diag([1, 0.5]); % Control cost
%optionally specify terminal cost matrix scaling, P_f = terminalRegionScaling*P_LQR (infinite-time LQR at final state)
terminalRegionScaling = 1; % Terminal constraint cost
                            %[TUNEABLE] increasing this would decrease the volume of the terminal matrix, P_f
                            %and hence increase the terminal cost term (improve tracking/convergence performance) 
                            % most probably, values greater than 1 would work

run("Step2_FeedbackControllerSynthesis.m");
disp('- - - - - - -'); disp(" ");

drawnow
%keyboard

%% Optionally, do Monte-Carlo rollouts to check whether the TVLQR is stabilizing

close all;
startTimeIndex = 1; %start time for the rollouts
startMaxPerturbation = 1; %a measure of max initial perturbations to state
                         %decrease this for a smaller initial set
%run("./utils/checkClosedLoop_MCRollouts.m");
%keyboard

%% [Optional] Load all the saved files for further analysis

clearvars; close all;
addpath('./lib/');

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');

%following functions can take in an additional (optional) argument
%specifiying the projection dimensions - for example: [1 2], [1 3], etc.
projectionDims_2D = [1 2];
plotOneLevelSet_2D(x_nom, P, projectionDims_2D); %axis auto or %axis normal

plotOneLevelSet_3D(x_nom, P);
daspect([1 1 2]);

%for k=1:length(time_instances)
%    disp(matrix_condition_number(P(:,:,k)));
%    disp(' ');
%end

<<<<<<< HEAD
keyboard
=======
%keyboard
>>>>>>> unicycle

%% Polynomialize system dynamics for SOS (algebraic) programming and compute dynamics of state-deviations (xbar)

run("Step3_getDeviationDynamics.m");
%Note: comment out lines 118-122 of the above script 
%      if you don't want to double-check that xbar_dot(0) = 0 at each t = t_k

disp('- - - - - - -'); disp(" ");

%% Time-conditioned invariant set analysis (with time-dependance)
<<<<<<< HEAD
disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
clearvars; close all; 
%keyboard;
=======

% disp('Computing time-sampled invariant set certificates using SOS programming..'); disp('Hit Continue or F5'); disp(' ');
% clearvars; close all; 
% 
% % Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% % Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
% %             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
% %             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
% 
% % hyperParameter sweep to find a good initial guess values for level set boundary value rho
% rho0Max = 0.2; rho0NumVals = 10;
% rho0Vals = linspace(1e-3, rho0Max, rho0NumVals); %rho0 = 0 doesn't make intuitive sense, hence a low value to begin with
% 
% cMax = 5; cNumVals = 2*cMax + 1;
% cVals = linspace(-cMax, cMax, cNumVals);
% 
% drawFlag = 0;
% usageMode = 'feasibilityCheck'; %run just for an initial feasibility check
% 
% %2D struct array to store results from the hyperparameter sweep
% resultStruct = struct('rho0', [], 'c', [], 'success', 0, 'inletVolume', [], 'outletVolume', []);
% results = repmat(resultStruct, rho0NumVals, cNumVals);
% 
% %save the initialisation parameters
% for i = 1:length(rho0Vals)
%     for j = 1:length(cVals)
%         results(i,j).rho0 = rho0Vals(i);
%         results(i,j).c = cVals(j);
%     end
% end
% 
% for i = 1:length(rho0Vals)
%     for j = 1:length(cVals)
%         clearvars -except drawFlag usageMode rho0Vals cVals i j results
%         close all;
% 
%         rhoInitialGuessConstant = rho0Vals(i);
%         rhoInitialGuessExpCoeff = cVals(j);
% 
%         try                               
%             run("Step4_computeTimeSampledInvarianceCertificates.m");
%         catch
%             disp('Could not find a successful initial guess to start the alternation scheme!');
%             results(i,j).success = 0;
% 
%             break
%         end
% 
%         %storing computed data of successful run for further analysis      
%         results(i,j).success = 1;
%         results(i,j).inletVolume = initialInletVolume;
%         results(i,j).outletVolume = initialOutletVolume;
% 
%         results(i,j)
%         %keyboard
%     end
% end

% save('./precomputedData/hyperparamaterSweep.mat', 'rho0Vals', 'cVals', 'results');

%% Run the SOS program (optimisation) once the feasibility check or hyperparameter sweep passes through

close all
drawFlag = 1; %change this to '0' if you don't want to plot results
>>>>>>> unicycle

%specify SOS program hyperparameters
maxIter = 1; %maximum number of iterations

<<<<<<< HEAD
% Exponentially evolving rho_guess: rhoGuess_k = rho_0 * exp(c*(t_k - tf)/(t0 - tf)) 
% Usage note: c > 0 --> exp decreasing rho_guess (shrinking funnel -- preferred)
%             c = 0 --> constant rho_guess       ("tube" -- somewhat ideal)
%             c < 0 --> exp increasing rho_guess (expanding funnel -- not-so ideal)
rhoInitialGuessConstant = 0.05; %[TUNEABLE] rho_0: decrease value if initial guess fails, 
                                % keep it greater than 0!
rhoInitialGuessExpCoeff = 2; %[TUNEABLE] c: increase value if initial guess fails

usageMode = 'feasibilityCheck'; %run just for an initial feasibility check
try                               
    run("Step4_computeTimeSampledInvarianceCertificates.m");
catch
    disp('Could not find a successful initial guess to start the alternation scheme!');
end

%optimize once the feasibility check passes through
close all
=======
%best initialisation parameters after performing the sweep
rhoInitialGuessConstant = 0.05; %[TUNEABLE] rho_0: decrease value if initial guess fails, 
                                % keep it greater than 0!
rhoInitialGuessExpCoeff = 1; %[TUNEABLE] c: decrease value if initial guess fails
%Note:
%for "longer" trajectories,  rho0 = 0.05 and c = 1 or 1.5 seem to give a good feasible initial guess
%for "shorter" trajectories, rho0 = 0.03 and c = 0 gives a feasible initial guess

>>>>>>> unicycle
usageMode = 'shapeOptimisation'; %will have to workshop a name for this!
run("Step4_computeTimeSampledInvarianceCertificates.m");

disp('- - - - - - -'); disp(" ");

%% [Optional] Plot computed funnels
<<<<<<< HEAD
tic
=======
%tic
>>>>>>> unicycle
close all

projectionDims_2D = [1 3]; projectionDims_3D = [1 2 3];
run("./utils/plottingScript.m");
<<<<<<< HEAD
toc
=======
%toc
>>>>>>> unicycle
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

% Define the system dynamics: unicycle
% x = [p_x; p_y; theta]
% u = [v; omega]
function f = unicycle_dynamics(x, u)

    % Unicycle dynamics: x_dot = f(x, u)
    f = [u(1) * cos(x(3)); ...  % x_dot
         u(1) * sin(x(3)); ...  % y_dot
         u(2)];                % theta_dot
end

% Computes the condition number of a matrix to check for numerical ill-conditioning:
% kappa < 1e2 is good (k=1 is perfect condition and k > 1e3 is ill-conditioned)
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
        error('Matrix must be square');
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