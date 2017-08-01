function [  x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = pointMass_dyngoal_dyn( x, u )
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

% Error dynamics
e_goal_dot  = x_goal_dot - x_rob_dot;
e_goal_dot_e = x_goal_dot_x + x_rob_dot_x;
e_goal_dot_u = x_goal_dot_u - x_rob_dot_u;

% Compound system
x_dot = [x_rob_dot ; e_goal_dot];
x_dot_x = [x_rob_dot_x zeros(4) ; zeros(4) e_goal_dot_e];
x_dot_u = [x_rob_dot_u ; e_goal_dot_u];


Gamma = zeros(length(x));
Gamma(5:8,5:8) = 1*eye(4);
Sigma = eye(length(x));

end

