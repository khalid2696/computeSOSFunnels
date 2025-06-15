classdef plottingFnsClass
    properties %variables
    end
    
    methods %functions
        function M_2d = project_ellipsoid_matrix_2D(obj, M, projection_dims)
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
        
        function M_3d = project_ellipsoid_matrix_3D(obj, M, projection_dims)
            % Input:
            % M: nxn matrix defining the n-dimensional ellipsoid x^T M x < 1
            % projection_dims: 3-element vector specifying which dimensions to project onto
            %                  (e.g., [1 2 3] for xyz-plane, [1 3 2] for xzy-plane)
            
            n = size(M, 1); %get the dimensionality of matrix M
        
            basisMatrix = zeros(n,3);
            basisMatrix(projection_dims(1),1) = 1;
            basisMatrix(projection_dims(2),2) = 1;
            basisMatrix(projection_dims(3),3) = 1;
        
            M_3d = inv(basisMatrix' *inv(M) * basisMatrix);
        end
        
        function plotEllipse(obj, center, ellipseMatrix)
            
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
        
        % Function to visualize a 3D ellipsoid
        function plotEllipsoid(obj, center, ellipsoidMatrix, color)
        
            if nargin < 4 %assume a default color
                color = 'blue';
            end
        
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
            surf(X_e, Y_e, Z_e, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        end
    end
end