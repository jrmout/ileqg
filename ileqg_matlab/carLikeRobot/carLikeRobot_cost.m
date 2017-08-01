function [ l, l_x, l_xx, l_u, l_uu, l_ux  ] = carLikeRobot_cost( x, u, t, x_tar)
% CARLIKEROBOT_COST 
% Cost penalizing quadratic measure w.r.t a desired final state

sx = length(x);
su = length(u);

Q = zeros(sx);


Q_final = 1000*eye(sx);
%Velocity
Q_final(3,3) = 0;
Q_final(4,4) = 0;
Q_final(5,5) = 0;

R = 1e-2 * eye(su);


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
