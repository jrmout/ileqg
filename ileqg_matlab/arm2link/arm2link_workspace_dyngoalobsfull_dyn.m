function [  x_dot, x_dot_x, x_dot_u, Gamma, Sigma] = arm2link_workspace_dyngoalobsfull_dyn( x, u )
% 2-link arm-dynamics
% 
% x = (x_rob1 x_rob2 x_rob1_dot x_rob2_dot e_goal1 e_goal2 e_goal1_dot egoa2_dot)'
% u = (tau1 tau2)'
% e_goal = (x_goal - x_rob)

x_rob = x(1:4);
e_goal = x(5:8);
e_obs = x(9:12);

% Compute next state and control derivative
[x_rob_dot, x_rob_dot_u] = simulate2linkArm(x_rob,u);


% Finite differences
eps = 1e-6;

x1_plus = x_rob + eps*[1 0 0 0]';
x1_minus = x_rob - eps*[1 0 0 0]';
x2_plus = x_rob + eps*[0 1 0 0]';
x2_minus = x_rob - eps*[0 1 0 0]';
x3_plus = x_rob + eps*[0 0 1 0]';
x3_minus = x_rob - eps*[0 0 1 0]';
x4_plus = x_rob + eps*[0 0 0 1]';
x4_minus = x_rob - eps*[0 0 0 1]';

[x1_plus_dot, ~] = simulate2linkArm(x1_plus,u);
[x1_minus_dot, ~] = simulate2linkArm(x1_minus,u);
[x2_plus_dot, ~] = simulate2linkArm(x2_plus,u);
[x2_minus_dot, ~] = simulate2linkArm(x2_minus,u);
[x3_plus_dot, ~] = simulate2linkArm(x3_plus,u);
[x3_minus_dot, ~] = simulate2linkArm(x3_minus,u);
[x4_plus_dot, ~] = simulate2linkArm(x4_plus,u);
[x4_minus_dot, ~] = simulate2linkArm(x4_minus,u);

e1 = x1_plus_dot - x1_minus_dot;
e2 = x2_plus_dot - x2_minus_dot;
e3 = x3_plus_dot - x3_minus_dot;
e4 = x4_plus_dot - x4_minus_dot;

x_rob_dot_x = (1/(2*eps))*[e1 e2 e3 e4];
x_rob_dot_x(1:2, :) = [zeros(2) eye(2)];


% Point-mass Goal dynamics
% Simulates a linear Dim-dimensional point-mass with no orientation
Dim = 2;
%% Mass and Damping matrix for the admittance
M = 1*eye(Dim);
D = 1*eye(Dim);
%% Continuous-time dynamics of a point mass
A_goal = [zeros(Dim) eye(Dim);zeros(Dim) -(M\D)];
B = [zeros(Dim); inv(M)];
A = A_goal;

% Goal dynamics
x_goal = e_goal + x_rob;
x_goal_dot = A_goal*x_goal;
x_goal_dot_x = A_goal;
x_goal_dot_u = zeros(4,2);

x_goal_dot_x = A_goal;
x_goal_dot_u = zeros(4,2);

% Goal error dynamics
e_goal_dot  = x_goal_dot - x_rob_dot;
e_goal_dot_x = x_goal_dot_x + x_rob_dot_x;
e_goal_dot_u = x_goal_dot_u - x_rob_dot_u;

% Obstacle dynamics
x_obs = e_obs + x_rob;
A_obs = A_goal;
x_obs_dot = A_obs*x_obs;
x_obs_dot_x = A_obs;
x_obs_dot_u = zeros(4,2);

% Obstacle error dynamics
e_obs_dot  = x_obs_dot - x_rob_dot;
e_obs_dot_x = x_obs_dot_x + x_rob_dot_x;
e_obs_dot_u = x_obs_dot_u - x_rob_dot_u;

% Compound system
x_dot = [x_rob_dot ; e_goal_dot ; e_obs_dot];
x_dot_x = [x_rob_dot_x zeros(4) zeros(4);...
           zeros(4) e_goal_dot_x  zeros(4);... 
           zeros(4)  zeros(4) e_obs_dot_x];
x_dot_u = [x_rob_dot_u ; e_goal_dot_u ; e_obs_dot_u];


Gamma1 = zeros(length(x));
Gamma2 = zeros(length(x));
% Noisy obstacle Dyn
Gamma1(9:12,9:12) = 0.3*eye(4);
Gamma2(5:8,5:8) = 0.1*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;

Sigma(:,:,1) = eye(12);
Sigma(:,:,2) = eye(12);

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

if (1-D^2) > 0 
    % If the point is in the manipulator's workspace
    % this assumes elbow down position
    theta(2) = atan2(-sqrt(1-D^2),D);
    theta(1) = atan2(x(2),x(1)) - atan2(l2*sin(theta(2)),l1+l2*cos(theta(2)));
else 
    % If the point is not in the manipulator's workspace
    theta(2) = Inf;
    theta(1) = Inf;
end

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