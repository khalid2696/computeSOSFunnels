clc; clearvars; close all

%% Example 1: Centered ellipsoids
M1 = [2, 0, 0; 0, 0.5, 0; 0, 0, 4];
x_c1 = [0; 0; 0];

M2 = [4, 0, 0; 0, 1, 0; 0, 0, 9];
x_c2 = [0; 0; 0];

% Example 2: Non-centered ellipsoids
M3 = [6, 0, 0; 0, 3, 0; 0, 0, 1];
x_c3 = [2; 1; 0];

M4 = [5, 1, 0; 1, 3, 0; 0, 0, 2];
x_c4 = [1; 2; -1];

% Example 3: Large ellipsoid containing a smaller one
M5 = [1, 0.1, 0; 0.1, 1, 0; 0, 0, 1];
x_c5 = [0; 0; 0];

M6 = [4, 0, 0; 0, 4, 0; 0, 0, 4];
x_c6 = [0.4; 0.1; -0.25];

% Example 4: Highly eccentric ellipsoid
M7 = [10, 0, 0; 0, 1, 0; 0, 0, 0.1];
x_c7 = [3; -1; 2];

M8 = [1, 0.2, 0; 0.2, 2, 0; 0, 0, 5];
x_c8 = [-1; 1; 1];

% Example 5: Random positive definite matrices
M9 = 0.3*[3, 0.5, 0.1; 0.5, 2, 0.3; 0.1, 0.3, 1.5];
x_c9 = [0; 1; -1];

M10 = [4, -0.1, 0; -0.1, 3, 0.2; 0, 0.2, 2];
x_c10 = [0; 0.5; -0.35];


load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/setInvarianceCertificates.mat')

funnelInlet = ellipsoidMatrices(:,:,1)/rhoScaling(1);
funnelOutlet = ellipsoidMatrices(:,:,end)/rhoScaling(end);

M5 = funnelInlet;
M6 = funnelOutlet;

x_c5 = x_nom(:,end) + [0.1 -0.1 0.3]';
x_c6 = x_nom(:,end);

%% Assigning matrix parameters
EllipsoidMatrix1 = M5; EllipsoidCenter1 = x_c5;
EllipsoidMatrix2 = M6; EllipsoidCenter2 = x_c6;

figure
plotEllipsoid(EllipsoidMatrix1, EllipsoidCenter1, 'blue');
plotEllipsoid(EllipsoidMatrix2, EllipsoidCenter2, 'red');

% SDP-check
tic
check_usingSDP = isCompossible_usingSDP(EllipsoidMatrix1, EllipsoidCenter1, EllipsoidMatrix2, EllipsoidCenter2)
toc

%% SOS Programming

tic
pvar x1 x2 x3;   % define pvar "indeterminates"
x = [x1 x2 x3]'; % define the state vector

options.solver = 'sedumi';


%initialise the program
prog = sosprogram(x);

outerEllipsoidCondition = 1 - (x-EllipsoidCenter1)'*EllipsoidMatrix1*(x-EllipsoidCenter1);
innerEllipsoidCondition = 1 - (x-EllipsoidCenter2)'*EllipsoidMatrix2*(x-EllipsoidCenter2);

[prog,sL] = sossosvar(prog,1); %scalar multiplier

%multiplierPolyDeg = 4; %polynomial multiplier [even] degree (if scalar is not sufficient)
%[prog,s1] = sospolyvar(prog,monomials(x,multiplierPolyDeg)); %polynomial multiplier 
%prog = sosineq(prog, s1); 

%inner ellipse condition implies the outer one
%so the polynomial largrange term gets multiplied to the inner condition
p =  (outerEllipsoidCondition) - sL*(innerEllipsoidCondition);
prog = sosineq(prog, p); 

[prog, sol_info] = sossolve(prog, options);

if ~((sol_info.dinf==1) || (sol_info.pinf==1) || sol_info.feasratio < 0)
    disp('Given ellipsoids are contained');
    s1_sol = sosgetsol(prog, sL) %should be positive
else
    disp('Given ellipsoids are not contained');
end
toc

%% Function definitions

% Function to check ellipsoid containment
% checks whether ellipsoid 2 (red) is within ellipsoid 1 (blue)
% or alternatively whether ellipsoid 1 (blue) contains ellipsoid 2 (red)
function check = isCompossible_usingSDP(M_1, xc_1, M_2, xc_2)
    
    % Computes matrices for ellipsoid 1 (F, g, h)
    [F_1, g_1, h_1] = generate_ellipsoid_params(M_1, xc_1);

    % Computes matrices for ellipsoid 2 (F, g, h)
    [F_2, g_2, h_2] = generate_ellipsoid_params(M_2, xc_2);

    % Construct the LHS and RHS matrices for the SDP
    LHS = [F_1, g_1; g_1', h_1];  % LHS matrix 
    RHS = [F_2, g_2; g_2', h_2];  % RHS matrix (scaled by lambda)
    
    % SDP setup using SeDuMi or Mosek
    % Variables: lambda > 0
    cvx_begin sdp quiet
        cvx_solver mosek     % Use SeDuMi, or replace with cvx_solver mosek for Mosek
        variable lambda(1) nonnegative;  % Decision variable (lambda > 0)
        % Matrix inequality: LHS <= lambda * RHS
        LHS <= lambda * RHS;
    cvx_end

    % Output the results
    if strcmp(cvx_status, 'Solved')
        %fprintf('Ellipsoid containment is satisfied with lambda = %.4f\n', lambda);
        check = 1;
    else
        %fprintf('Ellipsoid containment is NOT satisfied.\n');
        check = 0;
    end
end

% Function to visualize a 3D ellipsoid
function plotEllipsoid(M, x_c, color)
    % Check if M is positive definite
    if ~all(eig(M) > 0)
        error('Matrix M must be positive definite.');
    end

    % Generate grid points on a unit sphere
    [X, Y, Z] = sphere(50); % Sphere with 50 x 50 resolution
    P = [X(:), Y(:), Z(:)]'; % Points on the unit sphere (3 x N)

    % Transform the unit sphere into the ellipsoid
    % Ellipsoid equation: (x - x_c)' * M * (x - x_c) = 1
    % Transform: M^(-1/2) * P
    M_inv_sqrt = sqrtm(inv(M)); % Compute M^(-1/2)
    EllipsoidPoints = M_inv_sqrt * P; % Scale points

    % Translate the ellipsoid to the center x_c
    for i = 1:size(EllipsoidPoints, 2)
        EllipsoidPoints(:, i) = EllipsoidPoints(:, i) + x_c; % Add x_c to each column
    end

    % Reshape points for surface plot
    X_e = reshape(EllipsoidPoints(1, :), size(X));
    Y_e = reshape(EllipsoidPoints(2, :), size(Y));
    Z_e = reshape(EllipsoidPoints(3, :), size(Z));

    % Plot the ellipsoid
    surf(X_e, Y_e, Z_e, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
    hold on;
end


function [F, g, h] = generate_ellipsoid_params(M, x_c)
    % Generates the parameters F, g, and h for an ellipsoid representation
    % in the form x'Fx + 2g'x + h <= 0 given the matrix M and center x_c
    % from the center-form representation (x-x_c)'M(x-x_c) < 1
    
    F = M;
    g = -M*x_c;
    h = x_c'*M*x_c - 1;
end