function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = pointMass_cost( x, u, t, x_tar)
% 2D_2ndOrder_Quadratic target 
%   Detailed explanation goes here
sx = length(x);
su = length(u);

Q = 10000*eye(sx);


Q_final = 10000000*eye(sx);
%Velocity
Q_final(3,3) = 1;
Q_final(4,4) = 1;

R = 1e-4 * eye(su);


if (isnan(t))
    l = (x-x_tar)'*Q_final*(x-x_tar);
    l_x = 2*Q_final*(x-x_tar);
    l_xx = 2*Q_final;
else
    l = (x-x_tar)'*Q*(x-x_tar) + u'*R*u;
    l_x = 2*Q*(x-x_tar);
    l_xx = 2*Q;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

