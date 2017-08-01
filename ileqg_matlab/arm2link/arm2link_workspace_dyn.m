function [  x_dot, x_dot_x, x_dot_tau, Gamma, Sigma] = arm2link_workspace_dyn( x, tau )
% 2-link arm-dynamics
% 
%  x = (x1 x2 x1_dot x2_dot)'
% tau = (tau1 tau2)'

% Compute next state and control derivative
[x_dot, x_dot_tau] = simulate2linkArm(x,tau);


% Finite differences
eps = 1e-4;

x1_plus = x + eps*[1 0 0 0]';
x1_minus = x - eps*[1 0 0 0]';
x2_plus = x + eps*[0 1 0 0]';
x2_minus = x - eps*[0 1 0 0]';
x3_plus = x + eps*[0 0 1 0]';
x3_minus = x - eps*[0 0 1 0]';
x4_plus = x + eps*[0 0 0 1]';
x4_minus = x - eps*[0 0 0 1]';

[x1_plus_dot, ~] = simulate2linkArm(x1_plus,tau);
[x1_minus_dot, ~] = simulate2linkArm(x1_minus,tau);
[x2_plus_dot, ~] = simulate2linkArm(x2_plus,tau);
[x2_minus_dot, ~] = simulate2linkArm(x2_minus,tau);
[x3_plus_dot, ~] = simulate2linkArm(x3_plus,tau);
[x3_minus_dot, ~] = simulate2linkArm(x3_minus,tau);
[x4_plus_dot, ~] = simulate2linkArm(x4_plus,tau);
[x4_minus_dot, ~] = simulate2linkArm(x4_minus,tau);

e1 = x1_plus_dot - x1_minus_dot;
e2 = x2_plus_dot - x2_minus_dot;
e3 = x3_plus_dot - x3_minus_dot;
e4 = x4_plus_dot - x4_minus_dot;

x_dot_x = (1/(2*eps))*[e1 e2 e3 e4];
x_dot_x(1:2, :) = [zeros(2) eye(2)];
%x_dot_x(3:4, 1:2) = [zeros(2)];


Gamma = zeros(length(x));
Sigma = eye(length(x));

end

function [x_dot, x_dot_tau] = simulate2linkArm(x,tau)
% Model parameters
m1 = 1.4;
m2 = 1.1;
l1 = 0.3;
l2 = 0.33;
s1 = 0.11;
s2 = 0.16;
i1 = 0.025;
i2 = 0.045;
d1 = i1 + i2 + m2*l1^2;
d2 = m2 * l1 * s2;
d3 = i2;
b11 = 0.05;
b22 = 0.05;
b12 = 0.025;
b21 = 0.025;


theta = zeros(4,1);

D = (x(1)^2+x(2)^2-l1^2-l2^2)/(2*l1*l2);
    

% this assumes elbow down position
theta(2) = atan2(-sqrt(1-D^2),D);
theta(1) = atan2(x(2),x(1)) - atan2(l2*sin(theta(2)),l1+l2*cos(theta(2)));

Jac = [ (- l1* sin(theta(1)) - l2 * sin(theta(1) + theta(2)))   (-l2 * sin(theta(1) + theta(2))) ; ...
                (l1 * cos(theta(1)) + l2 * cos(theta(1) + theta(2)))  (l2 * cos(theta(1) + theta(2))) ];
            
Jac_inv = Jac^(-1);

Jac_inv_dot =[ theta(4)*(sin(theta(1) + theta(2))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1))) - (cos(theta(1) + theta(2))*(l1*cos(theta(1) + theta(2))*cos(theta(1)) + l1*sin(theta(1) + theta(2))*sin(theta(1))))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1)))^2) + (theta(3)*sin(theta(1) + theta(2)))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1))),                                                                        - theta(4)*(cos(theta(1) + theta(2))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1))) + (sin(theta(1) + theta(2))*(l1*cos(theta(1) + theta(2))*cos(theta(1)) + l1*sin(theta(1) + theta(2))*sin(theta(1))))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1)))^2) - (theta(3)*cos(theta(1) + theta(2)))/(l1*cos(theta(1) + theta(2))*sin(theta(1)) - l1*sin(theta(1) + theta(2))*cos(theta(1))) ;
theta(4)*(((l2*cos(theta(1) + theta(2)) + l1*cos(theta(1)))*(l1*l2*cos(theta(1) + theta(2))*cos(theta(1)) + l1*l2*sin(theta(1) + theta(2))*sin(theta(1))))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1)))^2 - (l2*sin(theta(1) + theta(2)))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1)))) - (theta(3)*(l2*sin(theta(1) + theta(2)) + l1*sin(theta(1))))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1))), theta(4)*(((l2*sin(theta(1) + theta(2)) + l1*sin(theta(1)))*(l1*l2*cos(theta(1) + theta(2))*cos(theta(1)) + l1*l2*sin(theta(1) + theta(2))*sin(theta(1))))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1)))^2 + (l2*cos(theta(1) + theta(2)))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1)))) + (theta(3)*(l2*cos(theta(1) + theta(2)) + l1*cos(theta(1))))/(l1*l2*cos(theta(1) + theta(2))*sin(theta(1)) - l1*l2*sin(theta(1) + theta(2))*cos(theta(1)))];


theta(3:4) =  Jac_inv*x(3:4);


% Lagrangian dynamics
% inertia 
% M = (d1 + 2*d2*cos(theta2)   d3 + d2*cos(theta2) ; d3 + d2*cos(theta2) d3);
% centripetal + coriolis
% C = d2*sin(theta2)*(- theta2_dot*(2*theta1_dot + theta2_dot) ; theta1_dot^2);
% joint friction
% B = (b11 b12 ; b21 b22);

% Dynamics -> x_dot = F(x) + G(x)u
% F(x) = (f1 f2 f3 f4)';
d_theta2 = d1*d3 - d3^2 - (d2*cos(theta(2)))^2;
f1 = theta(3);
f2 = theta(4);

f3_prev = -d2*d3*(theta(3) + theta(4))^2 * sin(theta(2)) ...
    - d2^2 * theta(3)^2 * sin(theta(2)) * cos(theta(2)) ...
    - d2 * (b21*theta(3) + b22*theta(4)) * cos(theta(2)) ...
    + (d3*b11 - d3*b21)*theta(3) + (d3*b12 - d3*b22)*theta(4);
f3 = (1/d_theta2)* f3_prev;

f4_prev = d2*d3*theta(4)*(2*theta(3) + theta(4))*sin(theta(2)) ...
    + d1*d2*theta(3)^2*sin(theta(2)) + d2^2 * (theta(3) + theta(4))^2 * sin(theta(2))*cos(theta(2)) ...
    + d2 * ((2*b21 -b11)*theta(3) + (2*b22 - b12)*theta(4)) * cos(theta(2)) ...
    + (d1*b21 - d3*b11)*theta(3) + (d1*b22 - d3*b12)*theta(4);
f4 = (1/d_theta2)*f4_prev;

G = (1/d_theta2) * [ 0 0 ; 0 0 ; d3 -(d3 + d2*cos(theta(2))) ; -(d3 + d2*cos(theta(2))) d1 + 2*d2*cos(theta(2)) ];

theta_dot = [f1 f2 f3 f4]' + G*tau;

x_dot = zeros(4,1);
x_dot(1:2) = Jac*theta_dot(1:2);

x_dot(3:4) = Jac*(theta_dot(3:4) - Jac_inv_dot*x_dot(1:2));

G(3:4, :) = Jac_inv'*G(3:4,:);
x_dot_tau = G;
end