function symbolic_taylor_expansion(order) %Order of Taylor expansion (should be >= 1)
    
    if nargin < 1 %by default third-order approximation for good function-fit
        order = 3; 
    end
    
    %dimensionality of state-space and input-space
    nx = 3; nu = 2;
    x = createSymbolicVector('x', nx); %state vector
    u = createSymbolicVector('u', nu); %input vector

    %pvar x1 x2 x3 u1 u2
    %x = [x1 x2 x3]'; u = [u1 u2]'; %column matrix by convention -- nx1
    
    % Define the system dynamics: x_dot = f(x, u)
    f = [
        u(1) * cos(x(3));   % x_dot
        u(1) * sin(x(3));   % y_dot
        u(2);               % theta_dot
    ];
 
    %define the variables on which the function depends
    vars = [x; u]; %f(x,u) := f(vars)

    % Define symbolic expansion point
    expansion_point_varSymbol = 'a';
    expansion_point_a_symbolic = createSymbolicVector(expansion_point_varSymbol, length(vars));
    
    %pvar a1 a2 a3 a4 a5
    %expansion_point_a_symbolic = [a1 a2 a3 a4 a5]'; %column matrix by convention -- nx1

    taylor_approx = taylor_expansion(f, vars, expansion_point_a_symbolic, order);
    % this expression will be in terms of vars and expansion point_a

    % Display the symbolic result
    disp('Third-order Taylor series approximation (symbolic expansion point):');
    disp(taylor_approx);

    % A sample expansion point with numeric values
    expansion_point_a_numeric = [0, 0, 0, 1, 0]; % Sample point: x=0, y=0, theta=0, v=1, omega=0
    
    taylor_approx_at_a = evaluate_expansion_at_point(taylor_approx, expansion_point_varSymbol, expansion_point_a_numeric);
    % this expression will be in terms of only vars

    % Display the result
    disp('Taylor series approximation at the specific expansion point:');
    disp(simplify(taylor_approx_at_a));
    
    %to get info on the variable type
    %whos taylor_approx taylor_approx_at_a
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

function taylor_approx_at_a = evaluate_expansion_at_point(taylor_approx, point_varSymbol, numeric_a)
    
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
