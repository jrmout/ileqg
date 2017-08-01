function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = pointMass_dyngoalobs_cost( x, u, t)
% 2D_2ndOrder_Quadratic target 
%   Detailed explanation goes here
sx = length(x);
su = length(u);


Q_goal = zeros(4);



R = 1e-4 * eye(su);



wobs_width = 1/0.005;


% Exponential functional to obstacle distance
x_rob = x(1:4); 
e_goal = x(5:8);
e_obs_pos = x(9:10);
e_obs_vel = x(11:12);


if isnan(t)
    wobs = 0.01;
    Q_goal(1,1) = 1;
    Q_goal(2,2) = 1;
    
    l_obs = wobs*exp(- 0.5 * e_obs_pos' * wobs_width * e_obs_pos);
    l12obs_x = - wobs_width * l_obs * e_obs_pos;
    l12obs_xx = wobs_width^2 * l_obs*(e_obs_pos * e_obs_pos') - l_obs * wobs_width * (eye(2));
    

    
    l_obs_x = [l12obs_x' 0 0]';
    l_obs_xx = [l12obs_xx zeros(2) ; zeros(2,4)];


    l_goal = e_goal'*Q_goal*e_goal;
    l_goal_x = 2*Q_goal*e_goal;
    l_goal_xx = 2*Q_goal;


    l = l_goal + l_obs;
    l_x = [zeros(4,1) ; l_goal_x ; l_obs_x];
    l_xx = [zeros(4,12); zeros(4) l_goal_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx];
    
    % Ensure cost's positive definiteness
    [V,D] = eig(l_xx);
    d = diag(D);
    d(d<0) = 0;
    l_xx = V*diag(d)*V';
else
    Q_goal(1,1) = 10;
    Q_goal(2,2) = 10;
    wobs = 0.1;
    l_obs = wobs*exp(- 0.5 * e_obs_pos' * wobs_width * e_obs_pos);
    l12obs_x = - wobs_width * l_obs * e_obs_pos;
    l12obs_xx = wobs_width^2 * l_obs*(e_obs_pos * e_obs_pos') - l_obs * wobs_width * (eye(2));
    
%     %Ensure cost's positive definiteness
%     [V,D] = eig(l12obs_xx);
%     d = diag(D);
%     d(d<0) = 0;
%     l12obs_xx = V*diag(d)*V';
    
    l_obs_x = [l12obs_x' 0 0]';
    l_obs_xx = [l12obs_xx zeros(2) ; zeros(2,4)];


    l_goal = e_goal'*Q_goal*e_goal;
    l_goal_x = 2*Q_goal*e_goal;
    l_goal_xx = 2*Q_goal;


    l = l_goal + l_obs;
    l_x = [zeros(4,1) ; l_goal_x ; l_obs_x];
    l_xx = [zeros(4,12); zeros(4) l_goal_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx];
    
    % Ensure cost's positive definiteness
    [V,D] = eig(l_xx);
    d = diag(D);
    d(d<0) = 0;
    l_xx = V*diag(d)*V';
    
    
    l_control = u'*R*u;
    
    l = l + l_control;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end

end

