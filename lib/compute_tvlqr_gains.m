function [K, P] = compute_tvlqr_gains(f, x, u, x_nom, u_nom, Q, R, Pf, dt, method)
%f: symbolic nonlinear dynamics function depending on symbolic x and symbolic u
    
    if nargin < 10
        method = 'continuous'; %by default, use continuous-time equations for better numerical results
    end

    N = size(x_nom, 2);
    nx = length(x);
    nu = length(u);

    P = zeros(nx, nx, N);
    K = zeros(nu, nx, N);

    [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u);

    if strcmp(method, 'discrete')
        % === Discrete-Time Riccati Recursion ===
        [Kf, ~] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom(:,end), u_nom(:,end), Q, R);
        K(:,:,N) = Kf;
        P(:,:,N) = Pf;

        for k = N-1:-1:1
            % Get the value of the nominal trajectory and input at t=k
            x0_k = x_nom(:,k);
            u0_k = u_nom(:,k);

            %Get the linearised matrices (Jacobians) at nominal trajectory
            A_k = double(subs(A_sym, [x; u], [x0_k; u0_k]));
            B_k = double(subs(B_sym, [x; u], [x0_k; u0_k]));

            % Backward Riccati recursion
            P_k = P(:,:,k+1);
            P(:,:,k) = P_k + dt*(Q + A_k'*P_k + P_k*A_k - P_k*B_k*((R \ B_k')*P_k));
            K(:,:,k) = (R \ B_k')*P(:,:,k);
            %the above matrices are actually P and K at t=k-1
            %but stored in the matrix at index 'k'

            % Previous implementation: might be incorrect!
            %discrete-time system matrices
            %Ad_k = eye(nx) + dt * A_k;
            %Bd_k = dt * B_k;
    
            % Backward Riccati recursion
            %P_k = P(:, :, k+1);
            %S = R + Bd_k' * P_k * Bd_k;
            %K(:, :, k) = S \ (Bd_k' * P_k * Ad_k); % K_k = inv(S) * (B' * P_k * A)
            %P(:, :, k) = Q + Ad_k' * P_k * Ad_k - Ad_k' * P_k * Bd_k * (S \ (Bd_k' * P_k * Ad_k));
        end

    elseif strcmp(method, 'continuous')
        % === Continuous-Time DRE (ODE-based) ===
        tf = (N-1)*dt;
        tspan = linspace(tf, 0, N);

        A_all = zeros(nx, nx, N);
        B_all = zeros(nx, nu, N);
        for k = 1:N
            A_all(:,:,k) = double(subs(A_sym, [x; u], [x_nom(:,k); u_nom(:,k)]));
            B_all(:,:,k) = double(subs(B_sym, [x; u], [x_nom(:,k); u_nom(:,k)]));
        end

        A_interp = @(t) interp1(tspan, permute(A_all, [3 1 2]), t, 'linear', 'extrap');
        B_interp = @(t) interp1(tspan, permute(B_all, [3 1 2]), t, 'linear', 'extrap');

        P0_vec = reshape(Pf, [], 1);
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [T_sol, P_vec_sol] = ode15s(@(t, Pvec) dre_rhs(t, Pvec, A_interp, B_interp, Q, R, nx), [tf 0], P0_vec, opts);

        % Interpolate solution back to N steps
        P_interp = interp1(T_sol, P_vec_sol, tspan, 'linear');
        for k = 1:N
            P_k = reshape(P_interp(k,:), nx, nx);
            P_k = 0.5 * (P_k + P_k'); % enforce symmetry
            P(:,:,k) = P_k;

            A_k = A_all(:,:,k);
            B_k = B_all(:,:,k);
            K(:,:,k) = (R \ B_k') * P_k;
        end

        % Flip to go forward in time: [t0 to tf]
        P = flip(P, 3);    % P is nxnxN
        K = flip(K, 3);    % K is mxnxN
    end
end

% Differential Riccati equation, returns Pdot in vector format
function P_dot_vec = dre_rhs(t, Pvec, A_interp, B_interp, Q, R, n)
    P = reshape(Pvec, n, n);
    A = squeeze(A_interp(t));
    B = squeeze(B_interp(t));

    P_dot = -(P*A + A'*P - P*B*(R\(B'))*P + Q);
    
    % Vectorize the nxn matrix to n^2 x 1 vector (to use it with in-built ode solvers) 
    P_dot_vec = reshape(P_dot, [], 1);
end

% Computes symbolic Jacobians A and B from the dynamics x_dot = f(x, u)
function [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u)
    A_sym = jacobian(f, x);
    B_sym = jacobian(f, u);
end

% An infinite horizon LQR to stabilize the closed loop system at the terminal state
function [K_f, S_f] = compute_lqr_gain_at_terminal_state(f, x, u, x_nom_f, u_nom_f, Q, R)
    [A_sym, B_sym] = compute_symbolic_jacobians(f, x, u);
    A_f = double(subs(A_sym, [x; u], [x_nom_f; u_nom_f]));
    B_f = double(subs(B_sym, [x; u], [x_nom_f; u_nom_f]));
    [K_f, S_f, ~] = lqr(A_f, B_f, Q, R);
end