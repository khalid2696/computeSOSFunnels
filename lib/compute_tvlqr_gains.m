
%f: symbolic nonlinear dynamics function depending on symbolic x and symbolic u
function [K, P] = compute_tvlqr_gains(f, x, u, x_nom, u_nom, Q, R, Qf, dt)
    % Computes TVLQR gains using backward Riccati recursion
    N = size(x_nom, 2); %getting the number of discrete time-steps
    nx = length(x);
    nu = length(u);
    
    P = zeros(nx, nx, N);
    K = zeros(nu, nx, N);
    
    [K_f, ~] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom(:,end), u_nom(:,end), Q, R);
    K(:,:,N) = K_f; %terminal gain
    P(:,:,N) = Qf;  %terminal cost matrix
    
    for k = N-1:-1:1
  
        [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u);
        
        % Get the value of the nominal trajectory and input at t=k
        x0_k = x_nom(:,k);
        u0_k = u_nom(:,k);

        %Get the linearised matrices (Jacobians) at nominal trajectory
        A_k = double(subs(A_sym, [x; u], [x0_k; u0_k]));
        B_k = double(subs(B_sym, [x; u], [x0_k; u0_k]));

        %discrete-time system matrices
        Ad_k = eye(nx) + dt * A_k;
        Bd_k = dt * B_k;

        % Backward Riccati recursion
        P_k = P(:, :, k+1);
        S = R + Bd_k' * P_k * Bd_k;
        K(:, :, k) = S \ (Bd_k' * P_k * Ad_k); % K_k = inv(S) * (B' * P_k * A)
        P(:, :, k) = Q + Ad_k' * P_k * Ad_k - Ad_k' * P_k * Bd_k * (S \ (Bd_k' * P_k * Ad_k));
    end
end

%an infinite horizon LQR to stabilize the closed loop system at the terminal state
function [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom_f, u_nom_f, Q, R)

    [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u);
    
    %Get the linearised matrices (Jacobians) at terminal state
    A_f = double(subs(A_sym, [x; u], [x_nom_f; u_nom_f]));
    B_f = double(subs(B_sym, [x; u], [x_nom_f; u_nom_f]));
    
    [K_f,S_f,~] = lqr(A_f,B_f,Q,R);
end

% Pre-computed linearized dynamics around the nominal state/inputs: x0,u0
% just for unicyle dynamics (useful for debugging)
% Usage: [A, B] = linearize_dynamics(x_nom(:, k), u_nom(:, k), dt);
% function [A, B] = linearize_dynamics(x0, u0, dt)
%     % Linearize unicycle dynamics and discretize
%     theta = x0(3);
%     A_ct = [0, 0, -u0(1) * sin(theta);
%             0, 0,  u0(1) * cos(theta);
%             0, 0,  0];
%     B_ct = [cos(theta), 0;
%             sin(theta), 0;
%             0,          1];
% 
%     A = eye(3) + dt * A_ct;
%     B = dt * B_ct;
% end

% Computes symbolic Jacobians A and B from the dynamics x_dot = f(x, u)
function [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u)
    
    % Continuous-time Jacobians
    A_sym = jacobian(f, x);
    B_sym = jacobian(f, u);
end