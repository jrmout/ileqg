function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = arm2link_state_obstacle_cost( x, u, t, x_tar, x_obs)
% arm2link_cost
% Arm model parameters
l1 = 0.3;
l2 = 0.33;

sx = length(x);
su = length(u);



R = 1e-3 * eye(su);


% EF task coordinates
taskCoord = [l1*cos(x(1)) + l2*cos(x(1) + x(2))   l1*sin(x(1)) + l2*sin(x(1) + x(2))]';

dtaskdx = [ (- l1* sin(x(1)) - l2 * sin(x(1) + x(2)))   (-l2 * sin(x(1) + x(2))) ; ...
                (l1 * cos(x(1)) + l2 * cos(x(1) + x(2)))  (l2 * cos(x(1) + x(2))) ];
    
dtaskdxdx = [ (- l1* cos(x(1)) - l2 * cos(x(1) + x(2)))   (-l2 * cos(x(1) + x(2))) ( -l2 * cos(x(1) + x(2)))   (- l2 * cos(x(1) + x(2))) ; ...
            (- l1 * sin(x(1)) - l2 * sin(x(1) + x(2)))  (- l2 * sin(x(1) + x(2))) ( - l2 * sin(x(1) + x(2)))  (- l2 * sin(x(1) + x(2))) ...
             ];


% Consider final cost
% cost_final = w*(taskCoord - x_tar)'*(taskCoord - x_tar) + exp(- 0.5 * (taskCoord - x_obs)' * Q_obs * (taskCoord - x_obs))

w = 500;
wobs_width = (0.001)^(-1);
wobs = 10;


if (isnan(t))
    % Quadratic function to goal x_tar
    e = (taskCoord - x_tar);
    l_goal = w*(e'*e);
            
    l12goal_x = 2 * w * (e' * (dtaskdx))';
    edxdx = e'*dtaskdxdx;
    l12goal_xx = 2 * w * ((dtaskdx)' * (dtaskdx) + [edxdx(1:2) ; edxdx(3:4)]);
    
    l_x_goal = [l12goal_x' 0 0]';
    l_xx_goal = [l12goal_xx zeros(2) ; zeros(2,4)];

    % Exponential functional to avoid obstacle
    e_obs = (x_obs - taskCoord);
    l_obs = wobs*exp(- 0.5 * e_obs' * wobs_width * e_obs);
    
    l12obs_x = wobs_width * l_obs*((dtaskdx)'* e_obs );
    eobsdxdx = e_obs'*dtaskdxdx;
    l12obs_xx = wobs_width * l_obs*((dtaskdx)' * e_obs )'*((dtaskdx)' * e_obs ) + l_obs * wobs_width * ((dtaskdx)' * -(dtaskdx) + [eobsdxdx(1:2) ; eobsdxdx(3:4)]');
    l_x_obs = [l12obs_x' 0 0]';
    l_xx_obs = [l12obs_xx zeros(2) ; zeros(2,4)];
        
    l = l_goal;
    l_x = l_x_goal + l_x_obs;
    l_xx = l_xx_goal + l_xx_obs;
else
    % Exponential functional to avoid obstacle
    e_obs = (x_obs - taskCoord);
    l_obs = wobs*exp(- 0.5 * e_obs' * wobs_width * e_obs);
    
    l12obs_x = wobs_width * l_obs*((dtaskdx)'* e_obs );
    eobsdxdx = e_obs'*dtaskdxdx;
    l12obs_xx = wobs_width * l_obs*((dtaskdx)' * e_obs )'*((dtaskdx)' * e_obs ) + l_obs * wobs_width * ((dtaskdx)' * -(dtaskdx) + [eobsdxdx(1:2) ; eobsdxdx(3:4)]');
    l_x_obs = [l12obs_x' 0 0]';
    l_xx_obs = [l12obs_xx zeros(2) ; zeros(2,4)];
    
    
    l =  u'*R*u + l_obs;
    l_x = l_x_obs;
    l_xx = l_xx_obs;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

