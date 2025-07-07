clc; clearvars; close all

%% Add directories
addpath('../lib/');

%% Load the precomputed files

load('../precomputedData/nominalTrajectory.mat');
load('../precomputedData/LQRGainsAndCostMatrices.mat');
load('../precomputedData/setInvarianceCertificates.mat')

%% Plotting the funnels in 2D

%following functions can take in an additional (optional) argument
%specifiying the projection dimensions - for example: [1 2], [1 3], etc.

if ~exist('projectionDims_2D','var')
    projectionDims_2D = [1 2];
end 

plotFunnel_2D(x_nom, ellipsoidMatrices, rhoScaling, projectionDims_2D);
plotInitialSet_2D(x_nom(:,1), inletRegion, projectionDims_2D);
plotFinalSet_2D(x_nom(:,end), outletRegion, projectionDims_2D);

%formatting
title('Invariant Ellipsoidal Sets along the nominal trajectory');
%axis normal %or alternatively daspect([1 0.5]);

%% Plot 3D funnel

%following functions can take in an additional (optional) argument
%specifiying the projection dimensions - for example: [1 2 3], [1 3 4], etc.

if ~exist('projectionDims_3D','var')
    projectionDims_3D = [1 2 3];
end

plotFunnel_3D(x_nom, ellipsoidMatrices, rhoScaling, projectionDims_3D);

funnelInlet = ellipsoidMatrices(:,:,1)/rhoScaling(1);
plotInitialSet_3D(x_nom(:,1), funnelInlet, projectionDims_3D);

funnelOutlet = ellipsoidMatrices(:,:,end)/rhoScaling(end);
plotFinalSet_3D(x_nom(:,end), funnelInlet, projectionDims_3D);

%formatting
title('3-D Invariant Ellipsoidal Sets');
%axis normal %or alternatively 
daspect([1 1 0.5]);

%% Function definitions

function plotFunnel_2D(x_nom, ellipsoidMatrix, rhoScaling, projectionDims)
    figure; hold on; grid on; axis equal;

    P = plottingFnsClass();
    
    if nargin < 4
        projectionDims = [1 2]; %if not specified, by default x-y projection
    end

    %plot ellipsoidal invariant sets in 2D
    for k=1:1:size(x_nom,2)
        M = ellipsoidMatrix(:,:,k)/rhoScaling(k);
        M_xy = P.project_ellipsoid_matrix_2D(M, projectionDims);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k)]';
        P.plotEllipse(center, M_xy);
    end 
    
    %nominal trajectory
    plot(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),'--b');

    %formatting
    %xlabel('p_x');    ylabel('p_y');
    xlabel(['x_{', num2str(projectionDims(1)), '}'])
    ylabel(['x_{', num2str(projectionDims(2)), '}'])
end

function plotInitialSet_2D(x_initial, inletRegion, projectionDims)
    
    P = plottingFnsClass();

    if nargin < 3
        projectionDims = [1 2]; %if not specified, by default x-y projection
    end

    M = inletRegion;
    M_xy = P.project_ellipsoid_matrix_2D(M, projectionDims);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_initial(projectionDims(1)) + rotated_ellipse(1, :), ...
         x_initial(projectionDims(2)) + rotated_ellipse(2, :), ...
         '-.', 'Color', 'green', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
end

function plotFinalSet_2D(x_final, outletRegion, projectionDims)
    
    P = plottingFnsClass();

    if nargin < 3
        projectionDims = [1 2]; %if not specified, by default x-y projection
    end
    
    M = outletRegion;
    M_xy = P.project_ellipsoid_matrix_2D(M, projectionDims);

    [eig_vec, eig_val] = eig(M_xy);
    
    theta = linspace(0, 2*pi, 100); % Parameterize ellipse
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(x_final(projectionDims(1)) + rotated_ellipse(1, :), ...
         x_final(projectionDims(2)) + rotated_ellipse(2, :), ...
         '-.', 'Color', 'red', 'LineWidth', 1.2) %'FaceAlpha', 0.3); for 'fill' function
end

function plotFunnel_3D(x_nom, ellipsoidMatrix, rhoScaling, projectionDims)
    figure; view(3);
    hold on; grid on; axis equal;

    P = plottingFnsClass();

    if nargin < 4
        projectionDims = [1 2 3]; %if not specified, by default x-y-z projection
    end
    
    %plot ellipsoidal invariant sets in 3D
    for k=1:1:size(x_nom,2)
        M = ellipsoidMatrix(:,:,k)/rhoScaling(k);
        M_xyz = P.project_ellipsoid_matrix_3D(M, projectionDims);
        center = [x_nom(projectionDims(1),k), x_nom(projectionDims(2),k), x_nom(projectionDims(3),k)]';
        P.plotEllipsoid(center, M_xyz);
    end 
    
    outlet = ellipsoidMatrix(:,:,end)/rhoScaling(end);
    outlet_xyz = P.project_ellipsoid_matrix_3D(outlet, projectionDims);
    center = [x_nom(projectionDims(1),end), x_nom(projectionDims(2),end), x_nom(projectionDims(3),end)]';
    P.plotEllipsoid(center, outlet_xyz, 'red');

    %nominal trajectory
    plot3(x_nom(projectionDims(1),:),x_nom(projectionDims(2),:),x_nom(projectionDims(3),:),'--b');
    
    %formatting
    title('Invariant Ellipsoidal Sets along the nominal trajectory');
    %xlabel('p_x');    ylabel('p_y');    zlabel('p_z');
    xlabel(['x_{', num2str(projectionDims(1)), '}'])
    ylabel(['x_{', num2str(projectionDims(2)), '}'])
    zlabel(['x_{', num2str(projectionDims(3)), '}'])
end

function plotInitialSet_3D(x_initial, inletRegion, projectionDims)

    P = plottingFnsClass();

    if nargin < 3
        projectionDims = [1 2 3]; %if not specified, by default x-y-z projection
    end

    inlet_xyz = P.project_ellipsoid_matrix_3D(inletRegion, projectionDims);
    center = [x_initial(projectionDims(1)), x_initial(projectionDims(2)), x_initial(projectionDims(3))]';
    P.plotEllipsoid(center, inlet_xyz, 'green');
end

function plotFinalSet_3D(x_final, outletRegion, projectionDims)

    P = plottingFnsClass();

    if nargin < 3
        projectionDims = [1 2 3]; %if not specified, by default x-y-z projection
    end

    outlet_xyz = P.project_ellipsoid_matrix_3D(outletRegion, projectionDims);
    center = [x_final(projectionDims(1)), x_final(projectionDims(2)), x_final(projectionDims(3))]';
    P.plotEllipsoid(center, outlet_xyz, 'red');
end