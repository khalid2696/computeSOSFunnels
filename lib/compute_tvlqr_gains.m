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
        end

    elseif strcmp(method, 'continuous')
        % === Continuous-Time DRE (ODE-based) ===
        t_fwd = (0:N-1)*dt;                       % forward grid, matches x_nom indexing
    
        A_fun = matlabFunction(A_sym, 'Vars', {x, u});
        B_fun = matlabFunction(B_sym, 'Vars', {x, u});
        A_all = zeros(nx, nx, N);
        B_all = zeros(nx, nu, N);
        for k = 1:N
            A_all(:,:,k) = A_fun(x_nom(:,k), u_nom(:,k));
            B_all(:,:,k) = B_fun(x_nom(:,k), u_nom(:,k));
        end
    
        % Interpolants on the FORWARD grid
        A_gi = griddedInterpolant(t_fwd(:), reshape(A_all, nx*nx, N).', 'linear');
        B_gi = griddedInterpolant(t_fwd(:), reshape(B_all, nx*nu, N).', 'linear');
        A_interp = @(t) reshape(A_gi(t), nx, nx);
        B_interp = @(t) reshape(B_gi(t), nx, nu);
    
        % Integrate backward from tf to 0, output at the grid points
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [~, P_vec_sol] = ode15s(@(t,Pv) dre_rhs(t, Pv, A_interp, B_interp, Q, R, nx), ...
                                fliplr(t_fwd), reshape(Pf,[],1), opts);
        P_vec_sol = flipud(P_vec_sol);            % row k now corresponds to t_fwd(k)
    
        P = zeros(nx, nx, N);
        K = zeros(nu, nx, N);
        for k = 1:N
            P_k = reshape(P_vec_sol(k,:), nx, nx);
            P_k = 0.5*(P_k + P_k.');
            P(:,:,k) = P_k;
            K(:,:,k) = R \ (B_all(:,:,k).' * P_k);
        end
    end
end

% Differential Riccati equation, returns Pdot in vector format
function P_dot_vec = dre_rhs(t, Pvec, A_interp, B_interp, Q, R, nx)
    P  = reshape(Pvec, nx, nx);
    A  = A_interp(t);  B = B_interp(t);
    P_dot_vec = -(P*A + A.'*P - P*B*(R\(B.'*P)) + Q);
    P_dot_vec = P_dot_vec(:);
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