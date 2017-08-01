function [ x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = carLikeRobot_dyngoal_dyn( x, u )
% CARLIKE_DYN 
% Dynamics of a simplified 2D car model

% x_rob = [x1, x2, theta, v, k]'
% x_goal = [x1, x2, theta, v, k]'
% x = [x_rob ; (x_goal - x_rob)];

% x1,x2 -> 2D position
% theta -> orientation
% v -> velocity
% k -> curvature

% u = [ v_dot k_dot ]
x_rob = x(1:5);
e_goal = x(6:10);

% Robot dynamics
x_rob_dot = [x(4)*cos(x(3)) x(4)*sin(x(3)) x(4)*x(5) u(1) u(2)]';
x_rob_dot_x = [0 0 -x(4)*sin(x(3)) cos(x(3)) 0 ;...
     0 0  x(4)*cos(x(3)) sin(x(3)) 0 ;...
     0 0  0              x(5)    x(4);...
     0 0  0              0         0 ;...
     0 0  0              0         0 ];
 
x_rob_dot_u = zeros(5,2);
x_rob_dot_u(4,1) = 1;
x_rob_dot_u(5,2) = 1;

% Goal dynamics
x_goal = e_goal + x_rob;
x_goal_dot = [x_goal(4)*cos(x_goal(3)) x_goal(4)*sin(x_goal(3)) x_goal(4)*x_goal(5) 0 0]';
A_goal = [0 0 -x_goal(4)*sin(x_goal(3)) cos(x_goal(3)) 0 ;...
     0 0  x_goal(4)*cos(x_goal(3)) sin(x_goal(3)) 0 ;...
     0 0  0              x_goal(5)    x_goal(4);...
     0 0  0              0         0 ;...
     0 0  0              0         0];
x_goal_dot_x = A_goal;
x_goal_dot_u = zeros(5,2);

% Goal error dynamics
e_goal_dot  = x_goal_dot - x_rob_dot;
e_goal_dot_x = x_goal_dot_x + x_rob_dot_x;
e_goal_dot_u = x_goal_dot_u - x_rob_dot_u;


% Compound system
x_dot = [x_rob_dot ; e_goal_dot];
x_dot_x = [x_rob_dot_x -x_rob_dot_x ; -x_rob_dot_x e_goal_dot_x];
x_dot_u = [x_rob_dot_u ; e_goal_dot_u];

Gamma = eye(10);
Sigma = eye(10);
end

