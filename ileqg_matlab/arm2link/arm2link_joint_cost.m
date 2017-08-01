function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = arm2link_joint_cost( x, u, t, x_tar)
% arm2link_cost
%   Detailed explanation goes here

sx = length(x);
su = length(u);

Q = zeros(sx);

Q_final = 10000*eye(sx);
Q_final(3,3) = 10;
Q_final(4,4) = 10;

R = 1 * eye(su);


% EF task coordinates
%taskCoord = [l1*cos(x(1)) + l2*cos(x(1) + x(2)) l1*sin(x(1)) + l2*sin(x(1) + x(2))];
%e = (x_tar - taskCoord);

%dl_dx1 = 2*e(1)*(l1 *(-sin(x(1))) + l2*(-sin(x(1) + x(2)) )) +


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

