function [  x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = pointMass_dyngoalobsfull_dyn( x, u )
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
x_goal_dot = A*x_goal + B*[cos(2*x_goal(1)) sin(2*x_goal(2))]';
A_goal = A + B*[-2*sin(2*x_goal(1)) 0 0 0; 0 2*cos(2*x_goal(2)) 0 0];

x_goal_dot_x = A_goal;
x_goal_dot_u = zeros(4,2);

% Goal error dynamics
e_goal_dot  = x_goal_dot - x_rob_dot;
e_goal_dot_x = x_goal_dot_x + x_rob_dot_x;
e_goal_dot_u = x_goal_dot_u - x_rob_dot_u;

% Obstacle dynamics
x_obs = e_obs + x_rob;
x_obs_dot = A*x_obs + B*[cos(10*x_obs(1)) sin(10*x_obs(2))]';
A_obs = A + B*[-10*sin(10*x_obs(1)) 0 0 0; 0 10*cos(10*x_obs(2)) 0 0];
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
       
% x_dot_x = [x_rob_dot_x zeros(4) zeros(4);...
%            zeros(4) e_goal_dot_x x_rob_dot_x;... 
%            zeros(4) x_rob_dot_x e_obs_dot_x];

x_dot_u = [x_rob_dot_u ; e_goal_dot_u ; e_obs_dot_u];


Gamma1 = zeros(length(x));
Gamma2 = zeros(length(x));
% Noisy plant dynamics
%Gamma(1:4,1:4) = eye(4);
% Noisy goal dynamics
%Gamma(5:8,5:8) = eye(4);
% Noisy obstacle dynamics
Gamma1(9:12,9:12) = 0.3*eye(4);
Gamma2(5:8,5:8) = 0.3*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;

Sigma(:,:,1) = eye(12);
Sigma(:,:,2) = eye(12);

end

