function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = arm2link_state_cost( x, u, t, x_tar)
% arm2link_cost
% Arm model parameters
l1 = 0.3;
l2 = 0.33;

sx = length(x);
su = length(u);

w = 10000;

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

if (isnan(t))
    e = (taskCoord - x_tar);
    l = w*(e'*e);
            
    l12_x = 2 * w * (e' * (dtaskdx))';
    
    edxdx = e'*dtaskdxdx;
    
    l12_xx = 2 * w * ((dtaskdx)' * (dtaskdx) + [edxdx(1:2) ; edxdx(3:4)]);
    
    l_x = [l12_x' 0 0]';
    l_xx = [l12_xx zeros(2) ; zeros(2,4)];
else
    
    l =  u'*R*u;
    l_x = zeros(sx,1);
    l_xx = zeros(sx,sx);

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

