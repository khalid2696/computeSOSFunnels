clc; clearvars; close all

%% Load the precomputed files

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/setInvarianceCertificates.mat')

%% Status message
%Quantities at our disposal now

% N                  : number of time samples                   : scalar N
% time_instances     : time horizon (sampled)                   : 1 x N
% x_nom              : nominal state trajectory                 : n_x x N
% u_nom              : feedforward input tape                   : n_u x N
% candidateV         : Lyapunov certificates of invariance      : 1 x N (cell)
% ellipsoidMatrices  : ellipsoids characterizing candidateV     : n_x x n_x x N 
% rhoScaling         : level-set boundary value                 : 1 x N
% terminalRegion     : Ellipsoidal goal region in BRS analysis  : n x n
% multiplierTerms    : polynomial multipliers from S-procedure  : 1 x N-1 (cell)

N = length(time_instances);
n = size(x_nom, 1); m = size(u_nom, 1); %state and input vector dimensionality

%% Parameters

tolerance = 1e-6;
%% Further analysis

%run("plottingScript.m");
%daspect([1 1 2]) %maintains aspect (scaling)

%% Ellipsoid sampling

k = randi(length(time_instances))

tempState = x_nom(:,k);
tempInvariantSet = ellipsoidMatrices(:,:,k)/rhoScaling(k);

numSamplePoints = 1000;
%add argument: 'surface' to the following function for sampling only on the ellipsoid surface
ellipsoidSamples = samplePointsFromEllipsoid(tempState, tempInvariantSet, numSamplePoints);

figure; hold on; grid on; view(3);
plotEllipsoid(tempState, tempInvariantSet);

maxVal = 0;
minVal = 1;

for i =1:size(ellipsoidSamples,2)
    thisPoint = ellipsoidSamples(:,i);
    plot3(thisPoint(1),thisPoint(2),thisPoint(3),'xy');
    
    lyapFnValue = (tempState - thisPoint)'*tempInvariantSet*(tempState - thisPoint);
    
    %for catching errors
    if lyapFnValue < 0 || lyapFnValue > 1 + tolerance
        lyapFnValue   
        error('Sampled state may not be in the interior')
    end

    %determining the max and min values for debugging purposes
    if lyapFnValue > maxVal
        maxVal = lyapFnValue;
    end

    if lyapFnValue < minVal
        minVal = lyapFnValue;
    end

end

disp(['Minimum Lyapunv function Value: ' num2str(minVal)]);
disp(['Maximum Lyapunv function Value: ' num2str(maxVal)]);
disp(' ');

%% Function definitions

% Sampling interior or surface of an Ellipsoid   
function samplePoints = samplePointsFromEllipsoid(center, ellipsoidMatrix, numSamplePoints, samplingRegion)
    
    if nargin < 4
        samplingRegion = 'interior';
    end
    
    % Eigen decomposition
    [Q, Lambda] = eig(ellipsoidMatrix);
    n = size(center,1);

    % Semi-axis lengths
    semi_axes_lengths = 1 ./ sqrt(diag(Lambda));
    
    % Generate random points on the n-dimensional unit sphere
    sphere_points = randn(n, numSamplePoints); % Random points
    sphere_points = sphere_points ./ vecnorm(sphere_points); % Normalize to lie on the unit sphere
    
    if strcmpi(samplingRegion,'interior')
        scaling = rand(n,numSamplePoints); %multiply by scaling between 0 and 1 for interior
    else
        scaling = ones(n,numSamplePoints); %multiply by scaling of 1 for surface
    end

    sphere_points = sphere_points .* scaling;

    % Transform points to the ellipsoid
    samplePoints = Q * diag(semi_axes_lengths) * sphere_points + center;
end

% Function to visualize a 3D ellipsoid
function plotEllipsoid(center, ellipsoidMatrix, color)

    if nargin < 3 %assume a default color
        color = 'blue';
    end

    % Check if M is positive definite
    %if ~all(eig(ellipsoidMatrix) > 0)
    %    error('Matrix M must be positive definite.');
    %end

    % Generate grid points on a unit sphere
    [X, Y, Z] = sphere(50); % Sphere with 50 x 50 resolution
    P = [X(:), Y(:), Z(:)]'; % Points on the unit sphere (3 x N)

    % Transform the unit sphere into the ellipsoid
    % Ellipsoid equation: (x - x_c)' * M * (x - x_c) = 1
    % Transform: M^(-1/2) * P
    M_inv_sqrt = sqrtm(inv(ellipsoidMatrix)); % Compute M^(-1/2)
    EllipsoidPoints = M_inv_sqrt * P; % Scale points

    % Translate the ellipsoid to the center x_c
    for i = 1:size(EllipsoidPoints, 2)
        EllipsoidPoints(:, i) = EllipsoidPoints(:, i) + center; % Add x_c to each column
    end

    % Reshape points for surface plot
    X_e = reshape(EllipsoidPoints(1, :), size(X));
    Y_e = reshape(EllipsoidPoints(2, :), size(Y));
    Z_e = reshape(EllipsoidPoints(3, :), size(Z));

    % Plot the ellipsoid
    surf(X_e, Y_e, Z_e, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    axis equal;
    xlabel('p_x'); ylabel('p_y'); zlabel('\theta');
    grid on;
    hold on;
end