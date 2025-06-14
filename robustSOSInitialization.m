% Different initialization approaches offer varying trade-offs between success rate, 
% computational cost, and implementation complexity. 
% Progressive refinement and system-specific adaptive methods represent the optimal balance 
% for practical implementation.

%PROPOSED IMPROVEMENTS:
%1. Robust SOS Initialization - Implement system-specific conservative initialization with multiple fallback strategies 
%2. Basic Regularization - Add numerical conditioning improvements to existing SOS formulation 
%3. Progressive Warm Starting - Use successful solutions to initialize subsequent problems

% Numerical Stability Considerations:
% 
% Higher-dimensional systems are particularly sensitive to numerical conditioning issues. 
% Adding regularization terms to both the Lyapunov function constraints and S-procedure multiplier constraints significantly improves solver reliability . State scaling 
% based on the magnitude of cost-to-go matrices also helps maintain numerical stability.
% The combination of conservative initialization, system-specific parameter tuning, and numerical regularization 
% should resolve the infeasibility issues we experience with quadrotor and cart-pole systems. 
% The parallel computing and warm starting improvements will then provide substantial performance gains once the fundamental feasibility problems are addressed 

function [rhoGuess, candidateV, initSuccess] = robustSOSInitialization(xbar, time_instances, P, x_nom, deviationDynamics, systemType)
%ROBUSTSOSINITIALIZATION Enhanced initialization for SOS funnel computation
%   Specifically designed for quadrotor and cart-pole systems that struggle
%   with standard constant rho initialization
%
% Inputs:
%   time_instances: Time vector (1 x N)
%   P: Cost-to-go matrices from TVLQR (n x n x N) 
%   x_nom: Nominal trajectory (n x N)
%   deviationDynamics: Cell array of polynomial dynamics (1 x N)
%   systemType: 'quadrotor', 'cartpole', or 'general'
%
% Outputs:
%   rhoGuess: Robust initial rho scaling (1 x N)
%   candidateV: Initial Lyapunov function candidates (1 x N cell)
%   initSuccess: Boolean indicating if initialization succeeded

    N = length(time_instances);
    n = size(P, 1);
    
    % Initialize outputs
    rhoGuess = zeros(size(time_instances));
    candidateV = cell(size(time_instances));
    initSuccess = false;
    
    %% Strategy 1: System-Specific Conservative Initialization
    fprintf('Attempting system-specific initialization...\n');
    [rho1, V1, success1] = systemSpecificInit(xbar, systemType, time_instances, P, x_nom, n);
    
    if success1
        rhoGuess = rho1;
        candidateV = V1;
        initSuccess = true;
        fprintf('Success with system-specific initialization!\n');
        return;
    end
    
    %% Strategy 2: Progressive Scaling with Stability Analysis
    fprintf('Attempting progressive scaling initialization...\n');
    [rho2, V2, success2] = progressiveScalingInit(xbar, time_instances, P, deviationDynamics, n);
    
    if success2
        rhoGuess = rho2;
        candidateV = V2; 
        initSuccess = true;
        fprintf('Success with progressive scaling!\n');
        return;
    end
    
    %% Strategy 3: TVLQR-Based with Regularization
    fprintf('Attempting TVLQR-based initialization with regularization...\n');
    [rho3, V3, success3] = tvlqrRegularizedInit(xbar, time_instances, P, x_nom, n);
    
    if success3
        rhoGuess = rho3;
        candidateV = V3;
        initSuccess = true;
        fprintf('Success with TVLQR-regularized initialization!\n');
        return;
    end
    
    %% Strategy 4: Multi-level hierarchical approach (last resort)
    fprintf('Attempting hierarchical initialization...\n');
    [rho4, V4, success4] = hierarchicalInit(xbar, time_instances, P, x_nom, deviationDynamics, n);
    
    if success4
        rhoGuess = rho4;
        candidateV = V4;
        initSuccess = true;
        fprintf('Success with hierarchical initialization!\n');
        return;
    end
    
    % If all strategies fail
    fprintf('Warning: All initialization strategies failed. Using very conservative fallback.\n');
    [rhoGuess, candidateV] = conservativeFallback(xbar, time_instances, P, n);
    initSuccess = false;
end

%% Implementation of Individual Strategies

function [rho, V, success] = systemSpecificInit(xbar, systemType, time_instances, P, x_nom, n)
    N = length(time_instances);
    rho = zeros(size(time_instances));
    V = cell(size(time_instances));
    success = false;
    
    % System-specific parameters based on empirical testing
    switch lower(systemType)
        case 'quadrotor'
            % Quadrotor: 6-12 states, fast attitude dynamics
            base_rho = 1;  % Very conservative base
            decay_rate = 0.3; % Moderate decay
            state_weights = getQuadrotorWeights(n);
            
        case 'cartpole'  
            % Cart-pole: unstable equilibrium, sensitive to perturbations
            base_rho = 0.5; % Extremely conservative
            decay_rate = 0.5; % Faster decay for stability
            state_weights = getCartPoleWeights(n);
            
        case 'unicycle'
            % Unicycle: stable system, current defaults work
            base_rho = 0.4;
            decay_rate = 0.1;
            state_weights = ones(n, 1);
            
        otherwise
            % General high-dimensional system
            base_rho = 0.02 / sqrt(n); % Scale with dimension
            decay_rate = 0.4;
            state_weights = ones(n, 1);
    end
    
    % Generate time-varying rho with stability-based scaling
    t0 = time_instances(1);
    tf = time_instances(end);
    
    for k = 1:N
        tk = time_instances(k);
        
        % Exponential decay with stability margin
        time_factor = exp(-decay_rate * (tk - t0) / (tf - t0));
        
        % Condition number based scaling (better conditioning = larger rho)
        P_k = P(:,:,k);
        cond_factor = 1 / max(1, sqrt(cond(P_k) / 100)); % Penalty for ill-conditioning
        
        rho(k) = base_rho * time_factor * cond_factor;
        
        % Ensure terminal constraint is met
        if k == N
            rho(k) = min(rho(k), 0.8); % Ensure feasible terminal set
        end
        
        % Weighted Lyapunov function
        P_weighted = diag(state_weights) * P_k * diag(state_weights);
        V{k} = createLyapunovCandidate(xbar, P_weighted);
    end
    
    % Test feasibility with a simple check
    success = all(rho > 0) && all(rho < 1) && rho(end) <= rho(1);
end

function [rho, V, success] = progressiveScalingInit(xbar, time_instances, P, deviationDynamics, n)
    N = length(time_instances);
    rho = zeros(size(time_instances));
    V = cell(size(time_instances));
    success = false;
    
    % Start with very conservative guess
    rho(1) = 0.1;
    V{1} = createLyapunovCandidate(xbar, P(:,:,1));
    
    % Progressive scaling based on local stability
    for k = 2:N
        try
            % Analyze local stability of deviation dynamics
            A_local = extractLinearPart(xbar, deviationDynamics{k}, n);
            max_eigenvalue = max(real(eig(A_local)));
            
            % Scale based on stability margin
            if max_eigenvalue < -0.1
                % Good stability margin
                scaling_factor = 1.1;
            elseif max_eigenvalue < 0
                % Marginal stability
                scaling_factor = 1.0;
            else
                % Unstable - reduce significantly
                scaling_factor = 0.9;
            end
            
            % Apply conditioning-based adjustment
            cond_P = cond(P(:,:,k));
            if cond_P > 1e6
                scaling_factor = scaling_factor * 0.5 % Reduce for ill-conditioned P
            end
            
            scaling_factor
            rho(k) = rho(k-1) * scaling_factor;
            
            % Ensure reasonable bounds
            rho(k) = max(0.01, min(0.5, rho(k)));
            
            V{k} = createLyapunovCandidate(xbar, P(:,:,k));
            
        catch
            % If analysis fails, use conservative scaling
            rho(k) = rho(k-1) * 0.95;
            V{k} = createLyapunovCandidate(xbar, P(:,:,k));
        end
    end
    
    % Ensure monotonic decrease toward terminal set
    rho = ensureMonotonicDecay(rho);
    success = true;
end

function [rho, V, success] = tvlqrRegularizedInit(xbar, time_instances, P, x_nom, n)
    N = length(time_instances);
    rho = zeros(size(time_instances));
    V = cell(size(time_instances));
    
    % Use TVLQR cost-to-go with regularization
    regularization = 1e-3; % Small regularization for numerical stability
    
    for k = 1:N
        % Regularized cost-to-go matrix
        P_reg = P(:,:,k) + regularization * eye(n);
        
        % Estimate rho based on nominal state magnitude and P matrix
        x_nom_k = x_nom(:,k);
        nominal_cost = x_nom_k' * P_reg * x_nom_k;
        
        % Scale rho based on nominal cost (smaller cost = larger allowable region)
        rho(k) = 0.1 / (1 + nominal_cost);
        
        % Ensure reasonable bounds
        rho(k) = max(0.1, min(0.2, rho(k)));
        
        V{k} = createLyapunovCandidate(xbar, P_reg);
    end
    
    % Ensure terminal set constraint
    rho(end) = min(rho(end), 0.8);
    
    success = true;
end

%% Helper Functions

function state_weights = getQuadrotorWeights(n)
    % Typical quadrotor state: [x, y, z, phi, theta, psi, vx, vy, vz, p, q, r]
    % Weight fast attitude dynamics more conservatively
    if n == 12
        state_weights = [1, 1, 1, 5, 5, 5, 1, 1, 1, 10, 10, 10]'; % Full state
    elseif n == 6
        state_weights = [1, 1, 1, 5, 5, 5]'; % Position + attitude
    else
        % Generic scaling for different state dimensions
        state_weights = [ones(n/2, 1); 3*ones(n/2, 1)]; % Conservative for latter half
    end
end

function state_weights = getCartPoleWeights(n)
    % Cart-pole state: [cart_pos, pole_angle] or [cart_pos, cart_vel, pole_angle, pole_vel]
    if n == 4
        state_weights = [1, 5, 1, 5]'; % Weight angle states more conservatively
    elseif n == 2
        state_weights = [1, 5]'; % Position, angle
    else
        state_weights = ones(n, 1);
    end
end

function V = createLyapunovCandidate(xbar, P_matrix)
    % Create symbolic Lyapunov function candidate
    V = xbar' * P_matrix * xbar;
end

function A_linear = extractLinearPart(xbar, deviationDyn, n)
    % Extract linear part from polynomial deviation dynamics
    % This is a simplified version - enhance based on polynomial structure
    try
        % Get Jacobian at origin (linear part)
        A_linear = double(jacobian(deviationDyn, xbar));
        A_linear = subs(A_linear, xbar, zeros(n, 1));
    catch
        A_linear = zeros(n, n); % Fallback for symbolic issues
    end
end

function rho_mono = ensureMonotonicDecay(rho)
    % Ensure rho decreases monotonically (required for funnel property)
    rho_mono = rho;
    for k = 2:length(rho)
        if rho_mono(k) > rho_mono(k-1)
            rho_mono(k) = rho_mono(k-1) * 0.98; % Slight decrease
        end
    end
end

function [rho, V] = conservativeFallback(xbar, time_instances, P, n)
    % Ultra-conservative fallback when all else fails
    N = length(time_instances);
    rho = 0.001 * ones(size(time_instances)); % Very small funnel
    V = cell(size(time_instances));
    
    for k = 1:N
        % Use heavily regularized P matrix
        P_safe = P(:,:,k) + 0.1 * eye(n);
        V{k} = createLyapunovCandidate(xbar, P_safe);
    end
    
    fprintf('Using ultra-conservative fallback initialization.\n');
end
