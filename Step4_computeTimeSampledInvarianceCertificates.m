%clc; clearvars; close all

warning('off','MATLAB:singularMatrix');

%% Loading the time instances, nominal trajectory, nominal input and feedback control gain

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');
load('./precomputedData/deviationDynamics.mat');

%% Status message
%Quantities at our disposal now

% N                  : number of time samples                   : scalar N
% time_instances     : time horizon (sampled)                   : 1 x N
% x_nom              : nominal state trajectory                 : n_x x N
% u_nom              : feedforward input tape                   : n_u x N
% K                  : Feedback gains (sampled)                 : n_u x n_x x N
% P                  : cost-to-goal matrix (sampled)            : n_x x n_x x N
% xbar               : pvar type state deviations               : n_x x 1 
% systemPolyDynamics : polynomialised system dynamics (sampled) : 1 x N cell 
% deviationDynamics  : state deviation dynamics (sampled)       : 1 x N cell

%Note: cost-to-go matrices, P will be used for candidate Lyapunov functions
%Warning: systemPolyDynamics may not be useful because it's in terms of syms type variables 

N = length(time_instances);
n = size(x_nom, 1); m = size(u_nom, 1); %state and input vector dimensionality

%% Hyper-parameters

% 1. Parameters pertaining to SOS Programming
options.solver = 'mosek'; %behind-the-scenes SDP solver of choice
                          %options -- mosek/sedumi
multiplierPolyDeg = 6; %polynomial multiplier of predefined degree
LyapunovFnDeg = 2;     %quadratic Lyapunov function
tolerance = 1e-6;
convergenceTolerance = 1e-2; %if less than 1 percent change

% 2. Parameters pertaining to Alternation Scheme (SCP) -- feasibility step and optimisation step             
if ~exist('maxIter','var')
    maxIter = 3; %maximum number of iterations
end
rhoStepUpValue = 0.02; % analagous to alpha in gradient descent
                       %[TUNEABLE] decrease this if you run into infeasibility 
                       %Default Value: 0.001

% 3. Parameters pertaining to defining terminal set/goal region
startScaling = 0.8; %[TUNEABLE] increase this for a larger initial region
                    %decrease if computed funnel is weirdly shaped
goalScaling = 0.8;  %Keep it less than 1


% 4. Parameters pertaining to initial guess of level-set boundary value, rho
% can include a binary search to determine a good initial guess -- between 0 and 1

if ~exist('rhoGuessChoice','var')
    rhoGuessChoice = 'const'; %options: 'const' or 'exp' (if const doesn't work out)
end 

% Option1: constant rho_guess -- rhoGuessChoice = 'const'
%[TUNEABLE] decrease value if initial guess fails, keep it less than 1!
if ~exist('rhoInitialGuessConstant','var')
    rhoInitialGuessConstant = 0.4;
end 

% Option2: Exponentially (quickly) increasing rho_guess -- rhoGuessChoice = 'exp'
% (only if constant doesn't work)
% rhoGuess_k = 1 * exp(-c*(t_k - tf)/(t0 - tf)), 
% t_k = tf --> rhoVal = 1; t_k = 0 --> rhoVal = some small number (1/e^c); 

%[TUNEABLE] increase value if initial guess fails                         
if ~exist('rhoInitialGuessExpCoeff','var')
    rhoInitialGuessExpCoeff = 1.5;
end
%% Get the scaling for initial guess of level set boundary value, rho

[rhoInitialGuess, candidateV] = getInitialRhoGuessAndCandidateV (time_instances, xbar, deviationDynamics, P, ...
                                                                    rhoInitialGuessConstant, rhoInitialGuessExpCoeff, ...
                                                                       rhoGuessChoice); 
                                                                    %options - 'const' (option1) or 'exp' (option2)
ellipsoidMatrices = P; %initial guess of ellipsoid matrices are the cost-to-go matrices from TVLQR

%% Define the intial and final regions

goalRegionEllipsoidMatrix = goalScaling*P(:,:,end)/rhoInitialGuess(end);

startRegionEllipsoidMatrix = (1/startScaling)*P(:,:,1)/rhoInitialGuess(1);
%startRegionEllipsoidMatrix = (1/startScaling)*goalRegionEllipsoidMatrix;
% to force the inlet to be larger in size than outlet (!!might not behave properly!!)

%% Plot guess funnel and start/end regions

%plot(time_instances, rhoInitialGuess);
%plotFunnel(x_nom, ellipsoidMatrices, ones(size(rhoInitialGuess)));

plotFunnel(x_nom, ellipsoidMatrices, rhoInitialGuess);
title('Initial guess funnel');
plotInitialSet(x_nom(:,1), startRegionEllipsoidMatrix);
plotFinalSet(x_nom(:,end), goalRegionEllipsoidMatrix);
%% The first feasibility check to see whether we're able to find polynomial Lagrange multipliers at all time instances (for our guessV and guessRho) 

[~, ~, infeasibilityStatus] = findPolynomialMultipliers(time_instances, xbar, deviationDynamics, candidateV, rhoInitialGuess, ...
                                                                        startRegionEllipsoidMatrix, goalRegionEllipsoidMatrix, ...
                                                                          multiplierPolyDeg, options, tolerance);

plotFunnel(x_nom, ellipsoidMatrices, rhoInitialGuess);
plotInitialSet(x_nom(:,1), startRegionEllipsoidMatrix);
plotFinalSet(x_nom(:,end), goalRegionEllipsoidMatrix);
disp(rhoInitialGuess');

if ~infeasibilityStatus
    title('Initial guess V and rho scaling (feasible)');
    disp('Feasibility check passed for given rho guess..');
else
    title('Initial guess V and rho scaling (infeasible)');
    error('Could not find a successful initial guess to start the alternation scheme!')
end

%keyboard

%% Alternation Loop -- proceed if we're able to get a feasible initial guess

currRhoScaling = rhoInitialGuess;
prevRhoScaling = NaN(size(time_instances)); %NaN values to start the iterations

for iter=1:maxIter
    
    clc
    % ------ L-Step: Finding polynomial Lagrange multipliers  ------ %
    % Feasibility check to see whether we're able to determine multiplier terms 
    % at all time instances for a given guess of level-set boundary (rho) values and candidate V 
    [~, multiplierTerms, infeasibilityStatus] = findPolynomialMultipliers(time_instances, xbar, deviationDynamics, candidateV, currRhoScaling, ...
                                                                            startRegionEllipsoidMatrix, goalRegionEllipsoidMatrix, ...
                                                                              multiplierPolyDeg, options, tolerance);
 
    if ~infeasibilityStatus
        disp('Feasibility-check step: passed'); disp(' ');
        
        %plotFunnel(x_nom, ellipsoidMatrices, currRhoScaling);
        %plotInitialSet(x_nom(:,1), startRegionEllipsoidMatrix);
        %plotFinalSet(x_nom(:,end), goalRegionEllipsoidMatrix);
        %title('A valid funnel certificate');
    end

    if infeasibilityStatus
        if iter == 1
            %plotFunnel(x_nom, P, currRhoScaling);
            %plotFinalSet(x_nom(:,1), P(:,:,1)/rho_Start);
            %title('Scaling of cost-go-matrices (unsuccessful)')
            error('Could not find a successful initial guess to start the alternation scheme!')
        else
            disp('Infeasibility! Exiting the alternation scheme..')
            disp(' ');
            iter
            break
        end
    end

    feasibleMultiplierTerms = multiplierTerms;
    %keyboard
    
    % ------ V-step: Finding Lyapunov function & level-set boundary value ------ %
    % Optimising the level-set volume through level-set boundary (rho) values and Lyapunov functions V
    % at all time instances, for the feasible multipliers obtained in the previous step 
    [~, sol_candidateVArray, sol_rhoValsArray, infeasibilityStatus] = ...
                                  findLyapFnAndLevelSetValues(time_instances, xbar, deviationDynamics, candidateV, multiplierTerms, ...
                                                                startRegionEllipsoidMatrix, goalRegionEllipsoidMatrix, ...
                                                                   LyapunovFnDeg, options, tolerance); 

    if infeasibilityStatus
        disp('Infeasibility! Exiting the alternation scheme..')
        disp(' ');
        iter
        break
    end
    
    % Extracting the solutions and converting to Matrix format
    ellipsoidMatrices = NaN(n,n,length(time_instances));
    currRhoScaling = NaN(size(time_instances));
    
    for k=1:length(time_instances)
        V_polyFn = sol_candidateVArray{k};
        ellipsoidMatrices(:,:,k) = getEllipsoidMatrix_nD(V_polyFn, n);
        currRhoScaling(k) = sol_rhoValsArray{k};

        disp(matrix_condition_number(ellipsoidMatrices(:,:,k)));
        disp(' ');
    end

    if ~infeasibilityStatus
        disp('Optimisation step: passed'); disp(' ');
        
        plotFunnel(x_nom, ellipsoidMatrices, currRhoScaling);
        plotInitialSet(x_nom(:,1), startRegionEllipsoidMatrix);
        plotFinalSet(x_nom(:,end), goalRegionEllipsoidMatrix);
        
        title('Optimised funnel certificate');
        
        % display some volumetric properties
        disp('Volume of computed inlet set: '); disp(1/sqrt(det((ellipsoidMatrices(:,:,1)/currRhoScaling(1))))); disp(' ');
        disp('Volume of input initial set: '); disp(1/sqrt(det(startRegionEllipsoidMatrix))); disp(' ');

        disp('Volume of input final set: '); disp(1/sqrt(det(goalRegionEllipsoidMatrix))); disp(' ');
        disp('Volume of computed outlet set: '); disp(1/sqrt(det((ellipsoidMatrices(:,:,end)/currRhoScaling(end))))); disp(' ');

    end

    %assigning values for next iteration
    prevRhoScaling = currRhoScaling;
    currRhoScaling = (1 + rhoStepUpValue)*prevRhoScaling;
    %currRhoScaling = rhoScaleIncrements.*prevRhoScalintocg;

    startRegionEllipsoidMatrix = startRegionEllipsoidMatrix/(1 + rhoStepUpValue); %increasing the inlet volume
    
    candidateV = sol_candidateVArray;

    drawnow;
    iter
    %keyboard;
    
end

%% Display some relevant metrics

disp('Volume ratio of computed inlet to input inlet set (should be > 1):')
disp(sqrt(det(startRegionEllipsoidMatrix))/sqrt(det((ellipsoidMatrices(:,:,1)/currRhoScaling(1)))));

disp('Volume ratio of computed outlet to input outlet set (should be < 1):')
disp(sqrt(det(goalRegionEllipsoidMatrix))/sqrt(det((ellipsoidMatrices(:,:,end)/currRhoScaling(end)))));

disp('Volume ratio of computed inlet to outlet:')
disp(sqrt(det((ellipsoidMatrices(:,:,end)/currRhoScaling(end))))/sqrt(det((ellipsoidMatrices(:,:,1)/currRhoScaling(1)))));

%% Status message

clc

disp('Finished computing ellipsoidal invariance-certificates around the nominal trajectory');
disp(' ');

%% Save the multipliers, caniddate V and level-set boundary value to a file

rhoScaling = prevRhoScaling;
multiplierTerms = feasibleMultiplierTerms;
inletRegion = startRegionEllipsoidMatrix;
outletRegion = goalRegionEllipsoidMatrix;
save('./precomputedData/setInvarianceCertificates.mat', 'time_instances', 'candidateV', 'rhoScaling', 'ellipsoidMatrices', 'inletRegion', 'outletRegion', 'multiplierTerms');

disp('Saved the time-sampled ellipsoidal matrices parametrizing the invariant sets to a file!');
disp(' ');
%% ---------------------  Function definitions  ----------------------------

%% SOS Program Functions -- the two alternating steps

function [prog, sol_multipliersArray, infeasibilityStatus] = ...
    findPolynomialMultipliers(time_instances, xbar, deviationDynamics, candidateV, rhoGuess, M_i, M_f, ...
                              multiplierPolyDeg, options, tolerance)
    
    N = length(time_instances);
    
    %N-1 decision polynomial terms: multipliers(k)
    lagrangeMultipliers = cell(N+1,1);  %this'll hold the SOS vars  (before solving) 
    sol_multipliersArray = cell(N+1,1); %this'll hold the solutions (after solving)

    infeasibilityStatus = 0;

    %initialise the SOS program
    prog = sosprogram(xbar);
    
    for k = 2:1:N
    %parfor k = 2:N
        
        %sampling time - Ts
        deltaT = time_instances(k) - time_instances(k-1);
        
        %calculate Vdot
        fbar = deviationDynamics{k};
        V = candidateV{k};
        partial_dVdt = (candidateV{k} - candidateV{k-1})/deltaT;
        Vdot = partial_dVdt + jacobian(V, xbar)*fbar; 
        
        %determine rho and rhoDot
        rho = rhoGuess(k);
        rhoDot = (rhoGuess(k) - rhoGuess(k-1))/deltaT;
        %rhoDot = rhoDotArray(k);

        %multiplier polynomial
        [prog, sL] = sospolyvar(prog,monomials(xbar,0:multiplierPolyDeg));
        lagrangeMultipliers{k} = sL;
        %or alternatively
        %[prog, s1] = sossosvar(prog,monomials(x,0:multiplierPolyDeg/2));
        %(USE ONLY FOR ROA! NOT FOR INVARIANT SET COMPUTATION)

        % SOS constraints

        % 1. Positive definiteness of V (taken care by construction of V)
        prog = sosineq(prog, V - tolerance*(xbar'*xbar)); 
    
        % 2. Vdot constraint (generalised S-procedure)
        prog = sosineq(prog, (rhoDot - Vdot) - sL*(rho - V) - tolerance*(xbar'*xbar));

        % 3. non-negativity of multiplier polynomial (USE ONLY FOR ROA! NOT FOR INVARIANT SET COMPUTATION)
        %prog = sosineq(prog, s1);
    end
    
    % 4. Inlet constraints -- inlet set contains the user-defined initial set
    outerEllipsoidCondition = rhoGuess(1) - candidateV{1};
    innerEllipsoidCondition = 1 - xbar'*M_i*xbar; %given initial (ellipsoid) set
    
    [prog,s0] = sossosvar(prog,1); %scalar multiplier
    lagrangeMultipliers{1} = s0;
    prog = sosineq(prog, outerEllipsoidCondition - s0*innerEllipsoidCondition);    
    
    % 5. Outlet constraints -- outlet set contained within the user-defined final set
    outerEllipsoidCondition = 1 - xbar'*M_f*xbar; %given final (ellipsoid) set
    innerEllipsoidCondition = rhoGuess(end) - candidateV{end}; 
    
    [prog,sN] = sossosvar(prog,1); %scalar multiplier
    lagrangeMultipliers{end} = sN;
    prog = sosineq(prog, outerEllipsoidCondition - sN*innerEllipsoidCondition - tolerance*(xbar'*xbar));  

    % Solve the feasibility program
    [prog, sol_info] = sossolve(prog, options);
    
    if (sol_info.dinf==1) || (sol_info.pinf==1) || sol_info.feasratio < 0
        %disp('Infeasible! Change rho guess value..');
        %disp(' '); 
        infeasibilityStatus = 1;
        return
    end

    %extract the solution and store it in a cell array
    for k= 2:1:N+1
        temp_sol = sosgetsol(prog,lagrangeMultipliers{k});
        %temp_sol = temp_sol/max(temp_sol.coefficient); %normalise the coefficients
        %temp_sol = cleanpoly(temp_sol, tolerance); %clean up the terms (remove coefficients smaller than tolerance)
        
        sol_multipliersArray{k} = temp_sol;
    end
end


function [prog, sol_candidateVArray, sol_rhoValsArray, infeasibilityStatus] = ...
    findLyapFnAndLevelSetValues(time_instances, xbar, deviationDynamics, Vguess, multiplierTerms, ...
                                M_i, M_f, LyapunovFnDeg, options, tolerance)
    
    N = length(time_instances);
    
    %N decision scalar variables: rho(k)
    rhoValArray = cell(N,1);      %this'll hold the SOS vars  (before solving)
    sol_rhoValsArray = cell(N,1); %this'll hold the solutions (after solving)

    %N decision scalar variables: rho(k)
    candidateVArray = cell(N,1);      %this'll hold the SOS vars  (before solving)
    sol_candidateVArray = cell(N,1);  %this'll hold the solutions (after solving)

    infeasibilityStatus = 0;

    %initialise the SOS program with indeterminates
    prog = sosprogram(xbar);
    
    %getting stuff ready
    for k = 1:N
        [prog, rho] = sossosvar(prog,1);
        rhoValArray{k} = rho;

        [prog, V] = sospolyvar(prog,monomials(xbar,2:LyapunovFnDeg));
        candidateVArray{k} = V;
    end
    
    %imposing SOS constraints
    objective = 0;
    objective = objective + rhoValArray{1}; %include rho corresponding to inlet
    for k = 2:1:N
        
        deltaT = time_instances(k) - time_instances(k-1);

        fbar = deviationDynamics{k};

        V = candidateVArray{k};
        partial_dVdt = (1/deltaT)*(candidateVArray{k} - candidateVArray{k-1}); 
        Vdot = partial_dVdt + jacobian(V, xbar)*fbar; 
        
        rho = rhoValArray{k};
        rhoDot = (1/deltaT)*(rhoValArray{k} - rhoValArray{k-1});
        %SOSTOOLS seems to work only with multiplication operator and not
        %division for 'dpvars', hence, 1/deltaT*(..) instead of (..)/deltaT !

        %SOS Constraints
        % 1. Positive definiteness of V (should be taken care by construction of V)
        prog = sosineq(prog, V - tolerance*(xbar'*xbar));
    
        % 2. Vdot constraint (generalised S-procedure)
        sL = multiplierTerms{k};
        prog = sosineq(prog, (rhoDot - Vdot) - sL*(rho - V) - tolerance*(xbar'*xbar));

        %Objective function
        objective = objective + rho;
    end

    % 3. Normalisation of coefficients of V (sum of coefficients is same as the guessed Lyapunov function (from prev step))
    for k = 1:N
        %prog = soseq(prog, subs(V,xbar,ones(size(xbar))) - 1);
        prog = soseq(prog, subs(candidateVArray{k},xbar,ones(size(xbar))) - subs(Vguess{k},xbar,ones(size(xbar))));
    end

    % 4. Inlet constraints -- inlet set contains the user-defined initial set
    outerEllipsoidCondition = rhoValArray{1} - candidateVArray{1};
    innerEllipsoidCondition = 1 - xbar'*M_i*xbar; %given initial (ellipsoid) set
    
    [prog,s0] = sossosvar(prog,1); %scalar multiplier
    prog = sosineq(prog, outerEllipsoidCondition - s0*innerEllipsoidCondition - tolerance*(xbar'*xbar));    
    
    % 5. Outlet constraints -- outlet set contained within the user-defined final set
    outerEllipsoidCondition = 1 - xbar'*M_f*xbar; %given final (ellipsoid) set
    innerEllipsoidCondition = rhoValArray{end} - candidateVArray{end}; 
 
    sN = multiplierTerms{end};
    prog = sosineq(prog, outerEllipsoidCondition - sN*innerEllipsoidCondition - tolerance*(xbar'*xbar));  
 
    %inner ellipse condition implies the outer one
    %so the polynomial largrange term gets multiplied to the inner condition 

    %alternative objective -- just minimise the outlet region
    %objective = rhoValArray{N}; 

    %Set the objective
    prog = sossetobj(prog, objective); % Objective: minimise sum of rho (hence positive sign)
    
    % Solve the SOS program
    [prog, sol_info] = sossolve(prog, options);
    
    if (sol_info.dinf==1) || (sol_info.pinf==1) || sol_info.feasratio < 0
        infeasibilityStatus = 1;
        %disp('Infeasible problem!');
        %disp(' ');
        return
    end

    % Extract the results from solved SOS Program and store it in a cell array
    for k= 1:1:N      
        temp_rhoSol = double(sosgetsol(prog, rhoValArray{k}));
        sol_rhoValsArray{k} = temp_rhoSol;
        
        temp_VSol = sosgetsol(prog, candidateVArray{k});
        sol_candidateVArray{k} = temp_VSol;
    end
    
    %scaling the coefficients of V and rho by the same factor
    %for k= 1:1:N
    %    sol_candidateVArray{k} = temp_VSol/sol_rhoValsArray{end};
    %    sol_rhoValsArray{k} = sol_rhoValsArray{k}/sol_rhoValsArray{end};
    %end
end


%% Helper functions

function [rhoGuess, candidateV] = getInitialRhoGuessAndCandidateV(time_instances, xbar, deviationDynamics, costToGoMatrices, rho_0, c, option)
    
    N = length(time_instances);
    t0 = time_instances(1); tf = time_instances(N);
    
    rhoGuess = NaN(size(time_instances));
    candidateV = cell(size(time_instances));
    
    % Defining a quadratic lyapunov function initial candidate, V = xbar^T P xbar
    % xbar = x - x_nom; %state deviations *off* the nominal trajectory
    for k = 1:N
        tk = time_instances(k);
        
        if (strcmpi(option, 'const'))
            rhoGuess(k) = rho_0;
        else
            rhoGuess(k) = rho_0 * exp(-c*(tk - tf)/(t0 - tf));
        end

        %getting some initial Lyapunov candidates
        f = deviationDynamics{k};
        A = jacobian(f,xbar); %symbolic A matrix
        A_at_origin = double(subs(A,xbar,zeros(size(xbar)))); %numeric A matrix
        
        if ~all(real(eig(A_at_origin)) < 0)
            disp(k); disp(eig(A_at_origin));
            disp('Nominal trajectory cannot be determined to be stable using indirect Lyapunov method');  
            error('Check the closed-loop system synthesis -- rework the nominal trajectory computation and TVLQR synthesis!')
        end
    
        candidateV{k} = xbar'*costToGoMatrices(:,:,k)*xbar;
        
    end

    %rhoGuess(end) = 1;

end


function M = getEllipsoidMatrix_nD(V_polyFn, n)
    
    M = NaN(n); %empty nxn matrix to hold the ellipsoid matrix
    
    if length(V_polyFn.coefficient) <  n*(n-1)/2
        return %if we're missing some coefficients the following algorithm can't work out
    end        % n(n-1)/2 is the number of elements in an upper triangle matrix

    k = 1;
    for i=1:n
        for j=i:n
            if i == j %diagonal terms
                M(i,j) = double(V_polyFn.coefficient(k));
            else      %off-diagonal terms
                M(i,j) = double(V_polyFn.coefficient(k))/2;
                M(j,i) = double(V_polyFn.coefficient(k))/2;
            end

            k = k+1;
        end
    end
    
    M = full(M);
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

%% Plotting functions
function plotFunnel(x_nom, ellipsoidMatrix, rhoScaling)
    figure
    hold on;
    grid on; 
    axis equal;
    
    for k=1:length(x_nom)
        M = ellipsoidMatrix(:,:,k)/rhoScaling(k);
        M_xy = project_ellipsoid_matrix(M, [1 2]);
        center = x_nom(:,k);
        plotEllipse(center, M_xy)
    end 
    
    title('Invariant Ellipsoidal Sets along the nominal trajectory');
    xlabel('p_x');
    ylabel('p_y');
    plot(x_nom(1,:),x_nom(2,:),'--b');
end

function plotInitialSet(x_initial, initialEllipsoid)

    M = initialEllipsoid;
    M_xy = project_ellipsoid_matrix(M, [1 2]);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_initial(1) + rotated_ellipse(1, :), ...
         x_initial(2) + rotated_ellipse(2, :), ...
         '-.g', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
end

function plotFinalSet(x_initial, finalEllipsoid)

    M = finalEllipsoid;
    M_xy = project_ellipsoid_matrix(M, [1 2]);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_initial(1) + rotated_ellipse(1, :), ...
         x_initial(2) + rotated_ellipse(2, :), ...
         '-.r', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
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
