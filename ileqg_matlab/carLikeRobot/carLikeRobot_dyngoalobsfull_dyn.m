function [ x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = carLikeRobot_dyngoalobsfull_dyn( x, u )
% CARLIKE_DYN 
% Dynamics of a simplified 2D car model

% x_rob = [x1, x2, theta, v, k]'
% x_goal = [x1, x2, theta, v, k]'
% x_obs = [x1obs, x2obs, x1obs_dot, x2obs_dot]'
% x = [x_rob ; (x_goal - x_rob); (x_obs - x_rob_taskspace)];

% x1,x2 -> 2D position
% theta -> orientation
% v -> velocity
% k -> curvature

% u = [ v_dot k_dot ]

x_rob = x(1:5);
e_goal = x(6:10);
e_obs = x(11:14);
x_rob_taskspace = [x(1) x(2) x(4)*cos(x(3)) x(4)*sin(x(3))]';

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

% Simulates a linear Dim-dimensional point-mass with no orientation
Dim = 2;

% Mass and Damping matrix for the point mass obstacle
M = 10*eye(Dim);
D = 1*eye(Dim);

% Continuous-time dynamics of a point mass
A = [zeros(Dim) eye(Dim);zeros(Dim) -(M\D)];

% Obstacle dynamics
x_obs = e_obs + x_rob_taskspace;
A_obs = A;
x_obs_dot = A_obs*x_obs;
x_obs_dot_x_obs = A_obs;
x_obs_dot_u = zeros(4,2);

% Dynamics of the carlike robot in euclidean space
x_rob_taskspace_dot = [x(4)*cos(x(3)) x(4)*sin(x(3)) -x(4)^2*sin(x(3))*x(5)+u(1)*cos(x(3)) x(4)^2*cos(x(3))*x(5)+u(1)*sin(x(3))]';
x_rob_taskspace_dot_x_rob_taskspace = [0 0 1 0 ; 0 0 0 1 ; 0 0 0 0; 0 0 0 0];

x_rob_taskspace_dot_x_rob = [0 0 -x(4)*sin(x(3)) cos(x(3)) 0 ; ...
                             0 0  x(4)*cos(x(3)) sin(x(3)) 0 ; ...
                             0 0 -x(4)^2*cos(x(3))*x(5)-u(1)*sin(x(3)) -2*x(4)*sin(x(3))*x(5) -x(4)^2*sin(x(3)) ; ...
                             0 0  -x(4)^2*sin(x(3))*x(5)+ u(1)*cos(x(3)) 2*x(4)*cos(x(3))*x(5) x(4)^2*cos(x(3))]; % 4x5


x_rob_dot_x_rob_taskspace = [0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 1 0 0 0]'; %5x4   Not necessary...
x_rob_taskspace_dot_u = zeros(4,2);
x_rob_taskspace_dot_u(3,1) = cos(x(3));
x_rob_taskspace_dot_u(4,2) = sin(x(3));

% Obstacle error dynamics
e_obs_dot  = x_obs_dot - x_rob_taskspace_dot;
e_obs_dot_x = x_obs_dot_x_obs + x_rob_taskspace_dot_x_rob_taskspace;
e_obs_dot_u = x_obs_dot_u - x_rob_taskspace_dot_u;


% Compound system
x_dot = [x_rob_dot ; e_goal_dot ; e_obs_dot];
% x_dot_x = [x_rob_dot_x zeros(5) zeros(5,4) ; ...
%            zeros(5) e_goal_dot_x x_rob_dot_x_rob_taskspace; ...
%            zeros(4,5) x_rob_taskspace_dot_x_rob e_obs_dot_x];
x_dot_x = [x_rob_dot_x zeros(5) zeros(5,4) ; ...
   zeros(5) e_goal_dot_x zeros(5,4); ...
   zeros(4,5) zeros(4,5) e_obs_dot_x];
x_dot_u = [x_rob_dot_u ; e_goal_dot_u ; e_obs_dot_u];


Gamma1 = zeros(length(x));
Gamma2 = zeros(length(x));
% Noisy obstacle Dyn
Gamma1(11:14,11:14) = 0.1*eye(4);
Gamma2(6:10,6:10) = 0.1*eye(5);
Gamma2(9,9) = 0;
Gamma2(10,10) = 0;
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;
Sigma(:,:,1) = eye(14);
Sigma(:,:,2) = eye(14);

end

