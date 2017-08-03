function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = pointMass_dyngoal_cost( x, u, t)
% 2D_2ndOrder_Quadratic target 
%   Detailed explanation goes here
sx = length(x);
su = length(u);

Q_final = zeros(sx);

Q_final(5,5) = 1;
Q_final(6,6) = 1;
Q_final(7,7) = 0.1;
Q_final(8,8) = 0.1;

Q = Q_final;

R = 1e-5 * eye(su);


if (isnan(t))
    l = x'*Q_final*x;
    l_x = 2*Q_final*x;
    l_xx = 2*Q_final;
else
    l = x'*Q*x + u'*R*u;
    l_x = 2*Q*x;
    l_xx = 2*Q;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

