function [  x_all_dot, x_all_dot_x, x_all_dot_tau, Gamma, Sigma] = arm2link_joint_dyngoal_dyn( x_all, tau )
% 2-link arm-dynamics
% 
% x_all = (theta1 theta2 theta1_dot theta2_dot x_goal1 x_goal2 x_goal_dot1 x_goal_dot2)'
% tau = (tau1 tau2)'

x = x_all(1:4);
x_goal = x_all(5:8);

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

% Lagrangian dynamics
% inertia 
% M = (d1 + 2*d2*cos(theta2)   d3 + d2*cos(theta2) ; d3 + d2*cos(theta2) d3);
% centripetal + coriolis
% C = d2*sin(theta2)*(- theta2_dot*(2*theta1_dot + theta2_dot) ; theta1_dot^2);
% joint friction
% B = (b11 b12 ; b21 b22);

% Dynamics -> x_dot = F(x) + G(x)u
% F(x) = (f1 f2 f3 f4)';
d_x2 = d1*d3 - d3^2 - (d2*cos(x(2)))^2;
f1 = x(3);
f2 = x(4);

f3_prev = -d2*d3*(x(3) + x(4))^2 * sin(x(2)) ...
    - d2^2 * x(3)^2 * sin(x(2)) * cos(x(2)) ...
    - d2 * (b21*x(3) + b22*x(4)) * cos(x(2)) ...
    + (d3*b11 - d3*b21)*x(3) + (d3*b12 - d3*b22)*x(4);
f3 = (1/d_x2)* f3_prev;

f4_prev = d2*d3*x(4)*(2*x(3) + x(4))*sin(x(2)) ...
    + d1*d2*x(3)^2*sin(x(2)) + d2^2 * (x(3) + x(4))^2 * sin(x(2))*cos(x(2)) ...
    + d2 * ((2*b21 -b11)*x(3) + (2*b22 - b12)*x(4)) * cos(x(2)) ...
    + (d1*b21 - d3*b11)*x(3) + (d1*b22 - d3*b12)*x(4);
f4 = (1/d_x2)*f4_prev;

G = (1/d_x2) * [ 0 0 ; 0 0 ; d3 -(d3 + d2*cos(x(2))) ; -(d3 + d2*cos(x(2))) d1 + 2*d2*cos(x(2)) ];


x_dot = [f1 f2 f3 f4]' + G*tau;


%%% x_dot_x
f3dx2 = -(1/(d_x2^2))*2*(d2*cos(x(2)))*d2*sin(x(2)) * f3_prev + (1/(d_x2)) * ( -d2*d3*(x(3) + x(4))^2 * cos(x(2)) ...
 - d2^2 * x(3)^2 * (cos(x(2))^2 - sin(x(2))^2 ) ...
 + d2 * (b21*x(3) + b22*x(4)) * sin(x(2)));

f4dx2 = -(1/(d_x2^2))*2*(d2*cos(x(2)))*d2*sin(x(2)) * f4_prev + (1/(d_x2)) * ( d2*d3*x(4)*(2*x(3) + x(4))*cos(x(2)) ...
 + d1*d2*x(3)^2*cos(x(2)) + d2^2 * (x(3) + x(4))^2 * (cos(x(2))^2 - sin(x(2))^2 ) ...
 - d2 * ((2*b21 -b11)*x(3) + (2*b22 - b12)*x(4)) * sin(x(2)));

f3dx3 = (1/d_x2)*( -d2*d3*2*(x(3) + x(4)) * sin(x(2)) ...
    - d2^2 * 2*x(3) * sin(x(2)) * cos(x(2)) ...
    - d2 * b21 * cos(x(2)) ...
    + (d3*b11 - d3*b21));

f4dx3 = (1/d_x2)*( d2*d3*x(4)*2*sin(x(2)) ...
    + d1*d2*2*x(3)*sin(x(2)) + d2^2 * 2*(x(3) + x(4)) * sin(x(2))*cos(x(2)) ...
    + d2 * (2*b21 -b11) * cos(x(2)) ...
    + (d1*b22 - d3*b12) );

f3dx4 = (1/d_x2)*( -d2*d3*2*(x(3) + x(4)) * sin(x(2)) ...
    - d2 * b22 * cos(x(2)) ...
    + (d3*b12 - d3*b22) );

f4dx4 = (1/d_x2)*( d2*d3*(2*x(3) + 2*x(4))*sin(x(2)) ...
    + d2^2 * 2*(x(3) + x(4)) * sin(x(2))*cos(x(2)) ...
    + d2 * (2*b22 - b12) * cos(x(2)) ...
    + (d1*b22 - d3*b12) );

g3dx2 = -(1/(d_x2^2))*2*(d2*cos(x(2)))*d2*sin(x(2)) * [d3 -(d3 + d2*cos(x(2)))] * tau + 1/(d_x2) * (tau(2)*d2*sin(x(2)));

g4dx2 = -(1/(d_x2^2))*2*(d2*cos(x(2)))*d2*sin(x(2)) * [ -(d3 + d2*cos(x(2))) d1 + 2*d2*cos(x(2)) ]*tau ...
    + 1/(d_x2) * (tau(1)*d2*sin(x(2)) - tau(2)*2*d2*sin(x(2)));

x_dot_x = [0 0 1 0 ; 0 0 0 1 ; 0 f3dx2 f3dx3 f3dx4 ; 0 f4dx2 f4dx3 f4dx4] ...
    + [0 0 0 0 ; 0 0 0 0 ; 0 g3dx2 0 0 ; 0 g4dx2 0 0 ];

%%% x_dot_tau
x_dot_tau = G;


%%% A point mass navigating in task space
Dim = 2;

% Mass and Damping matrix for the obstacle
M = 1*(eye(Dim));
D = 1*(eye(Dim));

% Discrete System Dynamics with sample time of TA (x_{k+1} = A*x_k + B*u)
A = [zeros(Dim) eye(Dim);zeros(Dim) -(M\D)];
x_goal_dot = A * x_goal;

x_goal_dot_x = A;

x_goal_dot_tau = zeros(4,2);

x_all_dot = [x_dot ; x_goal_dot];

x_all_dot_x = [x_dot_x zeros(4); zeros(4) x_goal_dot_x];

x_all_dot_tau = [x_dot_tau ; x_goal_dot_tau];


Gamma = zeros(length(x_all));
Sigma = eye(length(x_all));

end

