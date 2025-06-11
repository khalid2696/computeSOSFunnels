%clc; clearvars; close all

%% Add directories
addpath('./lib/');

%% Loading the time instances, nominal trajectory, nominal input and feedback control gain
load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');

%% Status message
%Quantities at our disposal now

% N              : number of time samples        : scalar N
% time_instances : time horizon (sampled)        : 1 x N
% x_nom          : nominal state trajectory      : n_x x N
% u_nom          : feedforward input tape        : n_u x N
% K              : Feedback gains (sampled)      : n_u x n_x x N
% P              : cost-to-goal matrix (sampled) : n_x x n_x x N
% dynamicsFnHandle : the function handle of the system dynamics
% f_sym          : system dynamics, f in terms of symbolic x and u: f(x,u) = xdot

%Note: cost-to-go matrices, P will be used for candidate Lyapunov functions

N = length(time_instances);
n = size(x_nom, 1); m = size(u_nom, 1); %state and input vector dimensionality

% A quadratic lyapunov function candidate, V = delta_x^T P delta_x
%delta_x = x - x_nom; %state deviations *off* the nominal trajectory
%candidateV = delta_x'*P*delta_x;

%total control input, u_k = u_nom - K_k delta_x
%u_k = u_nom(:, k) - K(:, :, k) * (x_traj(:, k) - x_nom(:, k));

%% Symbolic Taylor expansion for any general expansion point 'a'

%create symbolic variables of state and input vector based on the
%respective dimensionality
x = createSymbolicVector('x', n); %state vector
u = createSymbolicVector('u', m); %input vector
%x = [x1 x2 x3]'; u = [u1 u2]'; %column matrix by convention -- nx1

% Define the system dynamics: x_dot = f(x, u)
symbolicDynamics = [
    u(1) * cos(x(3));   % x_dot
    u(1) * sin(x(3));   % y_dot
    u(2);               % theta_dot
];

%define the variables on which the function depends
symVars = [x; u]; %f(x,u) := f(vars)

% Define symbolic expansion point
expansion_point_varSymbol = 'a';
expansion_point_a_symbolic = createSymbolicVector(expansion_point_varSymbol, length(symVars));
%expansion_point_a_symbolic = [a1 a2 a3 a4 a5]'; %column matrix by convention -- nx1

order = 3;
disp('Hang on..')
disp(['Polynomializing the system dynamics using Taylor expansion of order ' num2str(order)]);
disp(' ');

%taylor_approx = taylor_expansion(symbolicDynamics, symVars, expansion_point_a_symbolic, order);
% this expression will be in terms of vars and expansion point_a

%pre-computed polynomial-approximated dynamics (for quicker results)!!
%variable name: taylor_approx (symbolic - nx1 matrix)
load('./precomputedData/taylorApproxDynamicsSym.mat');

%taylor-expansion of uniycle dynamics, order = 3 (for debugging)!!
%lines: 72--80

% %re-assigning variables for ease of usage
% x1 = x(1); x2 = x(2); x3 = x(3); u1 = u(1); u2 = u(2);
% a1 = expansion_point_a_symbolic(1); a2 = expansion_point_a_symbolic(2); a3 = expansion_point_a_symbolic(3);
% a4 = expansion_point_a_symbolic(4); a5 = expansion_point_a_symbolic(5);
% 
% taylor_approx = ...
% [a4*cos(a3) - cos(a3)*(a4 - u1) + a4*sin(a3)*(a3 - x3) - 2*sin(a3)*(a4 - u1)*(a3 - x3) - (a4*cos(a3)*(a3 - x3)^2)/2 - (a4*sin(a3)*(a3 - x3)^3)/6 + (3*cos(a3)*(a4 - u1)*(a3 - x3)^2)/2;
%  a4*sin(a3) - sin(a3)*(a4 - u1) - a4*cos(a3)*(a3 - x3) + 2*cos(a3)*(a4 - u1)*(a3 - x3) + (a4*cos(a3)*(a3 - x3)^3)/6 - (a4*sin(a3)*(a3 - x3)^2)/2 + (3*sin(a3)*(a4 - u1)*(a3 - x3)^2)/2;
%                                                                                                                                                                                    u2];
%
%% Compute the polynomial-ized system dynamics at each nominal (state, input) pair

%for debugging purposes
%rand_timeInstance = randi(N);
%nom_state = x_nom(:, rand_timeInstance);
%nom_input = u_nom(:, rand_timeInstance);

disp('Computing polynomial system dynamics at each nominal state-input pair');
disp(' ');

systemPolyDynamics = cell(1,N);
parfor k = 1:N %length of time samples

    nom_state = x_nom(:, k);
    nom_input = u_nom(:, k);
    expansion_point_a_numeric = [nom_state' nom_input'];

    sym_polyDynamics_at_a = get_expansion_at_point(taylor_approx, 'a', expansion_point_a_numeric);
    sym_polyDynamics_at_a = vpa(sym_polyDynamics_at_a, 4); %cleanup upto 4 decimal places

    systemPolyDynamics{k} = sym_polyDynamics_at_a;
end

%% Compute the deviation dynamics in terms of pvar type deviation: xbar

% Converting expression from syms to pvar (to use SOSTOOLS functionalities)
%pvar x1 x2 x3
%xbar = [x1 x2 x3]'; %deviations --> nx1 vector
%xbar: pvar variables referring to state deviations *off* the nominal state 
xbar = createPvarVector('x', n);

disp('Computing state deviation dynamics at each nominal state-input pair');
disp(' ');

deviationDynamics = cell(1,N);
parfor k = 1:N %length of time samples
    
    nom_state = x_nom(:, k);
    nom_input = u_nom(:, k);
    controlGainMatrix = K(:, :, k);
    sym_polyDynamics_at_a = systemPolyDynamics{k};
    
    %deviation_dynamics: polynomial function in terms of xbar (deviations) ONLY
    %symVars --> (symbolic) [x; u]
    deviationDynamics{k} = computeDeviationDynamics(sym_polyDynamics_at_a, symVars, xbar, nom_state, nom_input, controlGainMatrix);
    
    %checking whether at zero deviation, deviation dynamics fn value is zero
    check = subs(deviationDynamics{k}, xbar, zeros(size(xbar)));
    if (norm(double(check)) ~= 0)
        disp('At zero deviation, deviation dynamics fn value is NOT zero. Check it!');
    end
end

%% Saving all the relevant data
save('./precomputedData/deviationDynamics.mat', 'time_instances', 'xbar', 'systemPolyDynamics', 'deviationDynamics');

disp('Saved the time-sampled dynamics of state-deviations to a file!');
disp(' ');

%% Fiunction defintions
function deviation_dynamics = computeDeviationDynamics(taylor_approx_at_a, symVars, xBar_pvar, x_nom, u_nom, controlGainMatrix)
    
    % computing xdot = f value at the nominal state and input values
    nomStateInputVals = [x_nom' u_nom'];
    fVal_at_nomStateInput = subs(taylor_approx_at_a,[symVars(1) symVars(2) symVars(3) symVars(4) symVars(5)], nomStateInputVals);
    fVal_at_nomStateInput = vpa(fVal_at_nomStateInput, 6); %cleanup upto 4 decimal places

    state = x_nom + xBar_pvar; %nominal state + deviation
    input = u_nom - controlGainMatrix*xBar_pvar;
    %note that u1 and u2 are not used since they depend on x1, x2, x3 as
    %[u1; u2] = controlGainMatrix*[x1 x2 x3]'

    pvarVars = [state', input'];  % Corresponding pvar variables
    
    pf = [];
    
    for i = 1:length(taylor_approx_at_a)
        [coeffs_f, monomials_f] = coeffs(taylor_approx_at_a(i), symVars);  % Extract terms
        
        % Reconstruct in pvar form
        pvar_poly = 0;
        for j = 1:length(coeffs_f)
            % Convert symbolic monomial to pvar-compatible monomial
            monomial_pvar = 1;
            for k = 1:length(symVars)
                exp_k = feval(symengine, 'degree', monomials_f(j), symVars(k));  % Get exponent
                monomial_pvar = monomial_pvar * pvarVars(k)^double(exp_k);
            end
            % Sum the terms with coefficients
            pvar_poly = pvar_poly + double(coeffs_f(j)) * monomial_pvar;
        end
        pf = [pf; pvar_poly];  % Store the result
    end
    
    % System dynamics in terms of state and input deviations (pvar variables)
    systemDynamics_pvar_exp = pf;

    %deviation dynamics
    deviation_dynamics = systemDynamics_pvar_exp - double(fVal_at_nomStateInput);

    %clean-up the small terms less than the tolerance
    tolerance = 1e-6;
    deviation_dynamics = cleanpoly(deviation_dynamics, tolerance);
end

function taylor_approx = taylor_expansion(f, vars, a, order)
    % Compute the Taylor expansion of multivariate function f up to the given order
    num_vars = length(vars);
    taylor_approx = sym(zeros(size(f))); % Initialize the Taylor expansion as a zero vector
    
    % Iterate over each term in the vector-valued function f
    for i = 1:length(f)

        % Initialize the Taylor expansion for the i-th component
        f_i_taylor = subs(f(i), vars, a); % Zeroth-order term (f(a))
        
        % Add higher-order terms iteratively
        for k = 1:order
            term = 0; % Initialize the term of order k
            
            % Iterate over all multi-index combinations for the current order
            multi_indices = nchoosek(repmat(1:num_vars, 1, k), k);
            for j = 1:size(multi_indices, 1)
                % Compute the partial derivative for the current multi-index
                partial_derivative = f(i);
                product_term = 1;
                for idx = multi_indices(j, :)
                    partial_derivative = diff(partial_derivative, vars(idx));
                    product_term = product_term * (vars(idx) - a(idx));
                end
                
                % Add the current term to the total
                term = term + subs(partial_derivative, vars, a) * product_term / factorial(k);
            end
            
            % Add the k-th order term to the Taylor expansion for f(i)
            f_i_taylor = f_i_taylor + term;
        end
        
        % Store the Taylor expansion for the i-th component
        taylor_approx(i) = simplify(f_i_taylor); % Simplify the expression
    end
end

function taylor_approx_at_a = get_expansion_at_point(taylor_approx, point_varSymbol, numeric_a)
    
    subs_map = struct(); %initialise an empty structure
    for i = 1:length(numeric_a)
        temp = [point_varSymbol, num2str(i)];
        subs_map.(temp) = numeric_a(i);
    end
    % Substitute the values into the Taylor series expansion
    taylor_approx_at_a = subs(taylor_approx, subs_map);
end



%inputs 
% varName: string -- name of vector
% n: int -- length of vector
function sym_vector = createSymbolicVector(vecSymbol, n)
    % Generate variable names dynamically
    varNames = arrayfun(@(i) [vecSymbol, num2str(i)], 1:n, 'UniformOutput', false);
    
    % Define these as symbolic variables
    syms(varNames{:});
    
    % Construct a symbolic vector from these variables
    sym_vector = sym(varNames, 'real');

    sym_vector = sym_vector'; %column matrix by convention -- nx1
end

%inputs 
% varName: string -- name of vector
% n: int -- length of vector
function pvar_vector = createPvarVector(vecSymbol, n)
    var_names = arrayfun(@(i) [vecSymbol num2str(i)], 1:n, 'UniformOutput', false);
    var_string = strjoin(var_names, ' ');
    
    % Create pvar in caller workspace
    evalin('caller', ['pvar ' var_string]);
    
    % Create the vector in caller workspace and return it
    evalin('caller', ['temp_pvar_vector = [' strjoin(var_names, '; ') '];']);
    pvar_vector = evalin('caller', 'temp_pvar_vector');
    evalin('caller', 'clear temp_pvar_vector'); % clean up
end