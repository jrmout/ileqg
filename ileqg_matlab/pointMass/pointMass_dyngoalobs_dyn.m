function [  x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = pointMass_dyngoalobs_dyn( x, u )
% Simulates a linear Dim-dimensional point-mass with no orientation
Dim = 2;

%% Mass and Damping matrix for the admittance
M = 1*eye(Dim);
D = 1*eye(Dim);

%% Continuous-time dynamics of a point mass
A = [zeros(Dim) eye(Dim);zeros(Dim) -(M\D)];
B = [zeros(Dim); inv(M)];

x_rob = x(1:4);
e_goal = x(5:8);
e_obs = x(9:12);

% Robot dynamics
x_rob_dot = A*x_rob + B*u;
x_rob_dot_x = A;
x_rob_dot_u = B;

% Goal dynamics
x_goal = e_goal + x_rob;
A_goal = A;
x_goal_dot = A_goal*x_goal;
x_goal_dot_x = A_goal;
x_goal_dot_u = zeros(4,2);

% Goal error dynamics
e_goal_dot  = x_goal_dot - x_rob_dot;
e_goal_dot_x = x_goal_dot_x + x_rob_dot_x;
e_goal_dot_u = x_goal_dot_u - x_rob_dot_u;

% Obstacle dynamics
x_obs = e_obs + x_rob;
A_obs = A;
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
           zeros(4) e_goal_dot_x zeros(4);... 
           zeros(4) zeros(4) e_obs_dot_x];
x_dot_u = [x_rob_dot_u ; e_goal_dot_u ; e_obs_dot_u];


Gamma = zeros(length(x));
% Noisy plant dynamics
%Gamma(1:4,1:4) = eye(4);
% Noisy goal dynamics
%Gamma(5:8,5:8) = eye(4);
% Noisy obstacle dynamics
Gamma(9:12,9:12) = eye(4);
Sigma = eye(length(x));
%Sigma(9:12,9:12) = eye(4);

end

