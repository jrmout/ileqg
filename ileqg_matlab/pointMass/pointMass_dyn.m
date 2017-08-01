function [  x_dot, x_dot_x, x_dot_u, Gamma, Sigma ] = pointMass_dyn( x, u )
% Simulates a linear Dim-dimensional point-mass with no orientation
Dim = 2;

%% Mass and Damping matrix for the admittance
M = 10*eye(Dim);
D = 1*eye(Dim);

%% Continuous-time dynamics of a point mass
A = [zeros(Dim) eye(Dim);zeros(Dim) -(M\D)];
B = [zeros(Dim); inv(M)];

x_dot = A*x + B*u;

x_dot_x = A;

x_dot_u = B;

Gamma = zeros(length(A));
Sigma = eye(length(A));

end

