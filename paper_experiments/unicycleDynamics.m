function xdot = unicycleDynamics(x, u)
% TRUE (non-polynomial) unicycle dynamics used by checkBoundaryValidity.m
% ASSUMED ORDER -- CHECK against your dynamics file:
%   state x = [x; y; theta],  input u = [v; omega]
xdot = [u(1)*cos(x(3));
        u(1)*sin(x(3));
        u(2)];
end
