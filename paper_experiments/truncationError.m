function [truncRel, errToMargin] = truncationError(runLog, x_nom, u_nom, Kgain, time_instances, ...
                                                   deviationDynamics, xbar, fTrue, Ms)
% truncRel(k)    = max over samples of Delta_k / |rhoDot_k|   (truncation error in Vdot, normalized)
% errToMargin(k) = max over samples of Delta_k / m_k           (< 1 => certificate holds under true dynamics)
% Delta_k = Vdot_true - Vdot_poly,  m_k = rhoDot_k - Vdot_poly >= 0
P = runLog.P; rho = runLog.rho; N = numel(time_instances); n = size(x_nom,1);
dT = time_instances(2) - time_instances(1); rng(0);
truncRel = nan(1,N); errToMargin = nan(1,N);
for k = 2:N
    Pk = P(:,:,k); Pkm = P(:,:,k-1); rhoDot = (rho(k) - rho(k-1))/dT;
    U = randn(n,Ms); U = U./vecnorm(U); Xb = sqrt(rho(k))*(sqrtm(Pk)\U);
    d = zeros(1,Ms); r = zeros(1,Ms);
    for s = 1:Ms
        xb    = Xb(:,s);
        fbar  = fTrue(x_nom(:,k)+xb, u_nom(:,k) - Kgain(:,:,k)*xb) - fTrue(x_nom(:,k), u_nom(:,k));
        fpoly = double(subs(deviationDynamics{k}, xbar, xb));
        delta = 2*xb'*Pk*(fbar - fpoly);
        mPoly = rhoDot - ((xb'*(Pk - Pkm)*xb)/dT + 2*xb'*Pk*fpoly);
        d(s) = delta/abs(rhoDot);  r(s) = delta/mPoly;
    end
    truncRel(k) = max(d); errToMargin(k) = max(r);
end
end