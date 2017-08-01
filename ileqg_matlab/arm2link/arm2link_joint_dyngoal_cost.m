function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = arm2link_state_dyngoal_cost( x_all, u, t)
% arm2link_cost
% Arm model parameters
l1 = 0.3;
l2 = 0.33;


sx = length(x_all);
su = length(u);

x = x_all(1:4);
x_goal = x_all(5:8);


x_tar = x_goal(1:2);

w = 100000;

R = 1 * eye(su);


% EF task coordinates
taskCoord = [l1*cos(x(1)) + l2*cos(x(1) + x(2))   l1*sin(x(1)) + l2*sin(x(1) + x(2))]';

dtaskdx = [ (- l1* sin(x(1)) - l2 * sin(x(1) + x(2)))   (-l2 * sin(x(1) + x(2))) ; ...
                (l1 * cos(x(1)) + l2 * cos(x(1) + x(2)))  (l2 * cos(x(1) + x(2))) ];
    
dtaskdxdx = [ (- l1* cos(x(1)) - l2 * cos(x(1) + x(2)))   (-l2 * cos(x(1) + x(2))) ( -l2 * cos(x(1) + x(2)))   (- l2 * cos(x(1) + x(2))) ; ...
            (- l1 * sin(x(1)) - l2 * sin(x(1) + x(2)))  (- l2 * sin(x(1) + x(2))) ( - l2 * sin(x(1) + x(2)))  (- l2 * sin(x(1) + x(2))) ...
             ];

% Consider final cost
% cost_final = w*(x_tar - taskCoord)'*(x_tar - taskCoord)

e = (taskCoord - x_tar);
l = w*(e'*e);
l12_x = 2 * w * (e' * (dtaskdx))';
edxdx = e'*dtaskdxdx;
l12_xx = 2 * w * ((dtaskdx)' * (dtaskdx) + [edxdx(1:2) ; edxdx(3:4)]);

l_x = zeros(8,1);
l_x(1:2) = l12_x;
l_xx = zeros(8);
l_xx(1:2, 1:2) = l12_xx;

if (~isnan(t))
    l =  l + u'*R*u;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

