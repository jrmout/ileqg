function [ x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = carLikeRobot_dyn( x, u )
% CARLIKE_DYN 
% Dynamics of a simplified 2D car model

% x = [x1, x2, theta, v, k]'

% x1,x2 -> 2D position
% theta -> orientation
% v -> velocity
% k -> curvature

% u = [ v_dot k_dot ]

x_dot = [x(4)*cos(x(3)) x(4)*sin(x(3)) x(4)*x(5) u(1) u(2)]';

x_dot_x = [0 0 -x(4)*sin(x(3)) cos(x(3)) 0 ;...
     0 0  x(4)*cos(x(3)) sin(x(3)) 0 ;...
     0 0  0              x(5)    x(4);...
     0 0  0              0         0 ;...
     0 0  0              0         0 ];
 
x_dot_u = zeros(5,2);
x_dot_u(4,1) = 1;
x_dot_u(5,2) = 1;

Gamma = eye(5);
Sigma = eye(5);

end

