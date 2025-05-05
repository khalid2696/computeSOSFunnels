clc; clearvars; close all

%% Define the system dynamics

load('./precomputedData/nominalTrajectory.mat');
load('./precomputedData/LQRGainsAndCostMatrices.mat');
load('./precomputedData/deviationDynamics.mat');

%t_k = randi(length(time_instances));
t_k = 25;
gammaGuess = 0.01; %reduce this if alternation is infeasible

f = deviationDynamics{t_k};

%% Lyapunov indirect method for stability analysis

%Linearise the system at origin
A = jacobian(f,xbar); %symbolic A matrix

A_at_origin = double(subs(A,xbar,zeros(size(xbar)))); %numeric A matrix

if all(real(eig(A_at_origin)) < 0)
    disp('Systems linearised matrix, A is Hurwitz, all poles on the RHP:')
    disp(eig(A_at_origin));
    disp('Origin is stable by indirect Lyapunov method');  
end

%% Lyapunov direct method for stability analysis

Q = eye(length(xbar));
%Q = diag([1 1 0.5]);

lyap_P = lyap(A_at_origin',Q); %solves the lyapunov equation, A^TP + PA = Q
%the transpose is required because MATLAB solves this eqn instead: A*P + P*A' + Q = 0

tvlqr_P = P(:,:, t_k);

candidate_P = lyap_P;
%% SOS Programming

%parameters
options.solver = 'sedumi'; %behind-the-scenes SDP solver of choice
                          %options -- mosek/sedumi/SDPT3
multiplierPolyDeg = 6; %polynomial multiplier of predefined degree
LyapunovFnDeg = 2;     %quadratic Lyapunov function
tolerance = 1e-6;
valueStepupRate = 0.07; %analagous to alpha in gradient descent
                        %decrease this if you run into infeasibility
convergenceTolerance = 0.01; %if less than 1 percent change

%requires a candidate V and initial guess of gamma
candidateV = xbar'*candidate_P*xbar; % A quadratic lyapunov function candidate, V = x^TPx

%%
[prog, s1, infeasibilityStatus] = findingMultipliersStep(xbar, f, candidateV, gammaGuess, multiplierPolyDeg, options, tolerance);
%keyboard

infeasibilityStatus
%% Display data for figuring out stuff

if ~infeasibilityStatus
    disp('A feasible invariant ellipsoidal set (x^T M x <= 1):');
    disp(candidate_P/gammaGuess)
end


disp(['Time instance, t = ' num2str(time_instances(t_k)) ' seconds']);
disp(' ');
%disp(x_nom(:,t_k)');
%disp(u_nom(:,t_k)');

%disp('Matrix from solving Lyapunov equation:');
%disp(lyap_P)

%disp('Cost-to-go matrix (a possible initial candidate from TVLQR analysis)');
%disp(tvlqr_P);


%[prog, V_sol, gamma_sol, ~] = findingLyapFnAndLevelSetValueStep(x, f, candidateV, s1, LyapunovFnDeg, options, tolerance);
%keyboard

%[prog, gamma_sol] = findingLevelSetValueStep(x, f, candidateV, s1, options, tolerance);
%keyboard

%% Function definitions

function [prog, s1_sol, infeasibilityStatus] = findingMultipliersStep(x, f, candidateV, gammaGuess, multiplierPolyDeg, options, tolerance)
    
    %Compute the time-derivative of candidate Lyapunov function
    Vdot = jacobian(candidateV,x)*f;

    s1_sol = NaN;
    infeasibilityStatus = 0;
        
    %initialise the SOS program
    prog = sosprogram(x);
    
    %multiplier polynomial
    %[prog, s1] = sossosvar(prog,monomials(x,0:multiplierPolyDeg/2));
    
    %or alternatively
    [prog,s1] = sospolyvar(prog,monomials(x,0:multiplierPolyDeg));
    %temp = s1;
    
    % 1. non-negativity of multiplier polynomial (USE ONLY FOR ROA! NOT FOR INVARIANT SET COMPUTATION)
    %prog = sosineq(prog, s1);

    % 2. Positive definiteness of V (taken care by construction of V)
    prog = sosineq(prog, candidateV - tolerance*(x'*x)); 
    
    % 3. Vdot constraint (generalised S-procedure)
    prog = sosineq(prog, -Vdot - s1*(gammaGuess - candidateV) - tolerance*(x'*x));
    
    %solve the program
    [prog, sol_info] = sossolve(prog, options);
    
    if (sol_info.dinf==1) || (sol_info.pinf==1) %|| sol_info.feasratio < 0
        disp('Infeasible! Reduce gamma guess value..');
        disp(' '); 
        infeasibilityStatus = 1;
        return
    end

    %get the solution
    %s1_sol = sosgetsol(prog,temp);
    s1_sol = sosgetsol(prog,s1);
    s1_sol = s1_sol/max(s1_sol.coefficient); %normalise the coefficients
    s1_sol = cleanpoly(s1_sol, tolerance); %clean up the terms (remove coefficients smaller than tolerance)
end

function [prog, gamma_sol] = findingLevelSetValueStep(x, f, candidateV, multiplierPoly, options, tolerance)
    
    %initialise the SOS program
    prog = sosprogram(x);
    [prog,gamma] = sossosvar(prog,1); %scalar decision variable -- level-set boundary value
    %[prog,candidateV] = sospolyvar(prog,monomials(x,2:LyapunovFnDeg)); 
    
    %Constraints
    % 1. Positive definiteness of V (should be taken care by construction of V)
    prog = sosineq(prog, candidateV - tolerance*(x'*x)); 
    
    % 2. Vdot constraint (generalised S-procedure)
    Vdot = jacobian(candidateV,x)*f;
    prog = sosineq(prog, -Vdot - multiplierPoly*(gamma - candidateV) - tolerance*(x'*x));
     
    %Objective function
    prog = sossetobj(prog, -gamma); % Objective: maximize gamma
    
    % Solve the SOS program
    [prog, sol_info] = sossolve(prog, options);
    
    if (sol_info.dinf==1) || (sol_info.pinf==1) || sol_info.feasratio < 0
        error('Infeasible problem');
    end
    
    % Extract the results from solved SOS Program
    gamma_sol = double(sosgetsol(prog, gamma));
end

function [prog, V_sol, gamma_sol, infeasibilityStatus] = findingLyapFnAndLevelSetValueStep(x, f, candidateV, multiplierPoly, LyapunovFnDeg, options, tolerance)
    
    infeasibilityStatus = 0;

    %initialise the SOS program
    prog = sosprogram(x);
    [prog,gamma] = sossosvar(prog,1); %scalar decision variable -- level-set boundary value
    [prog,V] = sospolyvar(prog,monomials(x,2:LyapunovFnDeg)); 
    
    %Constraints
    % 1. Positive definiteness of V (taken care by construction of V)
    prog = sosineq(prog, V - tolerance*(x'*x)); 
    
    % 2. Vdot constraint (generalised S-procedure)
    Vdot = jacobian(V,x)*f;
    prog = sosineq(prog, -Vdot - multiplierPoly*(gamma - V) - tolerance*(x'*x));
    
    % 3. Normalisation of coefficients of V (sum of coefficients is same as the guessed Lyapunov function (from prev step))
    prog = soseq(prog, subs(V,x,ones(size(x))) - 1);
    %prog = soseq(prog, subs(V,x,ones(size(x))) - subs(candidateV,x,ones(size(x))));
    
    %Objective function
    prog = sossetobj(prog, -gamma); % Objective: maximize gamma
    
    % Solve the SOS program
    [prog, sol_info] = sossolve(prog, options);
    
    if (sol_info.dinf==1) || (sol_info.pinf==1) || sol_info.feasratio < 0
        error('Infeasible problem'); disp(' ');
        infeasibilityStatus = 1;
        return
    end
    
    % Extract the results from solved SOS Program
    gamma_sol = double(sosgetsol(prog, gamma));
    V_sol = sosgetsol(prog, V);
end

%% Helper functions
function M = getEllipsoidMatrix(V_polyFn)
    M = [double(V_polyFn.coefficient(1)) double(V_polyFn.coefficient(2))/2; 
        double(V_polyFn.coefficient(2)/2) double(V_polyFn.coefficient(3))];
    
    M = full(M);
end

%% sample usage of SOSTOOLS toolbox

%prog = sosprogram(x,gamma); % Declare state variables in the SOS program
%usage: SOSP = sosprogram(vars, decisionVars)

%scalar Lagrange multiplier
%[prog,s1] = sossosvar(prog,1); %scalar multiplier


%[prog,s1] = sossosvar(prog,monomials(x,0:multiplierPolyDeg/2)); 
%NOTE: By default 'sossosvar' gives us an sos expression
%we don't need to impose sosineq(s1) or sosineq(s2)
% s(x) = z(x)^T Q z, z(x) --> vector of monomials, Q --> psd matrix (decision variables)
% that is why we need to give "multiplierPolyDeg/2" for degree of z(x) monomial vector

%if you don't want this (by default sos), you can instead use 'sospolvar'
%For instance,
%LyapunovFnDeg = 4;
%[prog,V] = sospolyvar(prog,monomials(x,2:LyapunovFnDeg)) 

%the following 2 are equivalent:
% -- [prog,gamma] = sossosvar(prog,1); %scalar
% -- dpvar gamma; prog = sosdecvar(prog, gamma);

%% List of SOSTOOLS functions

% Monomial vectors construction:
%    MONOMIALS   --- Construct a vector of monomials with prespecified 
%                    degrees.
%    MPMONOMIALS --- Construct a vector of multipartite monomials with
%                    prespecified degrees.
%
% General purpose sum of squares program (SOSP) solver:
% SOSPROGRAM        --- Initialize a new SOSP.
% SOSDECVAR         --- Declare new decision variables in an SOSP.
% SOSPOLYVAR        --- Declare a new polynomial variable in an SOSP.
% SOSSOSVAR         --- Declare a new sum of squares variable in an SOSP.
% SOSPOLYMATRIXVAR  --- Declare a new matrix of polynomial variables in an SOSP.
% SOSSOSMATRIXVAR   --- Declare a new matrix of sum of squares polynomial
%                       variables in an SOSP.
% SOSPOSMATR        --- Declare a new positive semidefinite matrix variable in an SOSP
% SOSPOSMATRVAR     --- Declare a new symbolic positive semidefinite matrix 
%                       variable in an SOSP
% SOSQUADVAR        --- Declare a polynomial/SOS decision variable with
%                       customized structure.
% SOSEQ             --- Add a new equality constraint to an SOSP.
% SOSINEQ           --- Add a new inequality constraint to an SOSP.
% SOSMATRIXINEQ     --- Add a new matrix inequality constraint to an SOSP.
% SOSSETOBJ         --- Set the objective function of an SOSP.
% SOSPSIMPLIFY      --- Perform simplification of an SOSP.
% SOSSOLVE          --- Solve an SOSP.
% SOSGETSOL         --- Get the solution from a solved SOSP.
%
% Customized functions:
%    FINDSOS     --- Find a sum of squares decomposition of a given polynomial.
%    FINDLYAP    --- Find a Lyapunov function for a dynamical system.
%    FINDBOUND   --- Find a global/constrained lower bound for a polynomial.
