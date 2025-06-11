%clc; clearvars; close all

%% Add directories
addpath('../lib/');

%% Load the precomputed files

load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/LQRGainsAndCostMatrices.mat');
load('../precomputedData/setInvarianceCertificates.mat')

%% Plotting the projected funnels

plotFunnel(x_nom, ellipsoidMatrices, rhoScaling);
plotInitialSet(x_nom(:,1), inletRegion);
plotFinalSet(x_nom(:,end), outletRegion);

%% Function definitions

function plotFunnel(x_nom, ellipsoidMatrix, rhoScaling)
    figure; hold on; grid on; %axis equal;
    
    projection_dims = [1 3]; %x-theta space

    for k=1:1:length(x_nom)
        M = ellipsoidMatrix(:,:,k)/rhoScaling(k);
        M_xy = project_ellipsoid_matrix(M, projection_dims);
        center = [x_nom(projection_dims(1),k), x_nom(projection_dims(2),k)]';
        plotEllipse(center, M_xy)
    end 
    
    title('Invariant Ellipsoidal Sets along the nominal trajectory');
    xlabel('p_x');
    ylabel('\theta');
    plot(x_nom(projection_dims(1),:),x_nom(projection_dims(2),:),'--b');
end

function plotInitialSet(x_initial, initialEllipsoid)
    
    M = initialEllipsoid; 
    projection_dims = [1 3]; %x-theta space

    M_xy = project_ellipsoid_matrix(M, projection_dims);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_initial(projection_dims(1)) + rotated_ellipse(1, :), ...
         x_initial(projection_dims(2)) + rotated_ellipse(2, :), ...
         '-.g', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
end

function plotFinalSet(x_final, finalEllipsoid)

    M = finalEllipsoid;
    projection_dims = [1 3]; %x-theta space

    M_xy = project_ellipsoid_matrix(M, projection_dims);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_final(projection_dims(1)) + rotated_ellipse(1, :), ...
         x_final(projection_dims(2)) + rotated_ellipse(2, :), ...
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
    ellipseCenter = center; % 2D center of the ellipsoid
    [eig_vec, eig_val] = eig(ellipseMatrix);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(ellipseCenter(1) + rotated_ellipse(1, :), ...
         ellipseCenter(2) + rotated_ellipse(2, :), ...
         '-k', 'LineWidth', 1.2);  
end

% % Function to visualize a 3D ellipsoid
% function plotEllipsoid(center, ellipsoidMatrix, color)
% 
%     if nargin < 3 %assume a default color
%         color = 'blue';
%     end
% 
%     % Check if M is positive definite
%     %if ~all(eig(ellipsoidMatrix) > 0)
%     %    error('Matrix M must be positive definite.');
%     %end
% 
%     % Generate grid points on a unit sphere
%     [X, Y, Z] = sphere(50); % Sphere with 50 x 50 resolution
%     P = [X(:), Y(:), Z(:)]'; % Points on the unit sphere (3 x N)
% 
%     % Transform the unit sphere into the ellipsoid
%     % Ellipsoid equation: (x - x_c)' * M * (x - x_c) = 1
%     % Transform: M^(-1/2) * P
%     M_inv_sqrt = sqrtm(inv(ellipsoidMatrix)); % Compute M^(-1/2)
%     EllipsoidPoints = M_inv_sqrt * P; % Scale points
% 
%     % Translate the ellipsoid to the center x_c
%     for i = 1:size(EllipsoidPoints, 2)
%         EllipsoidPoints(:, i) = EllipsoidPoints(:, i) + center; % Add x_c to each column
%     end
% 
%     % Reshape points for surface plot
%     X_e = reshape(EllipsoidPoints(1, :), size(X));
%     Y_e = reshape(EllipsoidPoints(2, :), size(Y));
%     Z_e = reshape(EllipsoidPoints(3, :), size(Z));
% 
%     % Plot the ellipsoid
%     surf(X_e, Y_e, Z_e, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
%     axis equal;
%     xlabel('p_x'); ylabel('p_y'); zlabel('\theta');
%     grid on;
%     hold on;
% end