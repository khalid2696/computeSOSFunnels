function [viol, maxExcess, violPoly] = checkBoundaryValidity(runLog, x_nom, u_nom, Kgain, time_instances, ...
                                                            deviationDynamics, xbar, fTrue, Ms, MsPoly)
% Evaluates the sampled invariance condition (paper: eq. Lyapunov_invariance_condition_discretized)
% on points of the funnel boundary {V_k = rho_k}, k = 2..N:
%   (V_k - V_{k-1})/dT + grad V_k' * fbar  <=  (rho_k - rho_{k-1})/dT
% viol(k)      : fraction of Ms boundary samples violating it with the TRUE deviation dynamics
% maxExcess(k) : largest violation (Vdot - rhoDot) at knot k (<= 0 means no violation found)
% violPoly(k)  : same fraction with the Taylor polynomial (sanity check, expect ~0)
% Directions are sampled uniformly on the unit sphere and mapped to the ellipsoid boundary
% (uniform in direction, NOT uniform in surface area).

P   = runLog.P;  rho = runLog.rho;
N   = numel(time_instances);  n = size(x_nom, 1);
dT  = time_instances(2) - time_instances(1);
rng(0);

viol = nan(1, N);  maxExcess = nan(1, N);  violPoly = nan(1, N);
for k = 2:N
    Pk = P(:,:,k);  Pkm = P(:,:,k-1);
    rk = rho(k);    rhoDot = (rk - rho(k-1))/dT;

    U  = randn(n, Ms);  U = U./vecnorm(U);
    Xb = sqrt(rk)*(sqrtm(Pk)\U);                         % xb' * Pk * xb = rho_k

    ex = zeros(1, Ms);
    for s = 1:Ms
        xb   = Xb(:, s);
        fbar = fTrue(x_nom(:,k) + xb, u_nom(:,k) - Kgain(:,:,k)*xb) - fTrue(x_nom(:,k), u_nom(:,k));
        ex(s) = (xb'*(Pk - Pkm)*xb)/dT + 2*xb'*Pk*fbar - rhoDot;
    end
    viol(k)      = mean(ex > 0);
    maxExcess(k) = max(ex);

    if MsPoly > 0
        exP = zeros(1, MsPoly);
        for s = 1:MsPoly
            xb    = Xb(:, s);
            fpoly = double(subs(deviationDynamics{k}, xbar, xb));
            exP(s) = (xb'*(Pk - Pkm)*xb)/dT + 2*xb'*Pk*fpoly - rhoDot;
        end
        violPoly(k) = mean(exP > 0);
    end
end
end
