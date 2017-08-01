function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = pointMass_dyngoalobsfull_cost( x, u, t, sobs_pos, sobs_Sigma)
% 2D_2ndOrder_Quadratic target 
%   Detailed explanation goes here
sx = length(x);
su = length(u);


Q_goal = zeros(4);
R = 1e-4 * eye(su);
wobs_width = 1/0.002;


% Exponential functional to obstacle distance
x_rob = x(1:4); 
e_goal = x(5:8);
e_obs_pos = x(9:10);

% Static obstacles
Q_sobs = 1*eye(2);
n_sobs = size(sobs_pos, 2);

if isnan(t)
    wobs = 0.001;
    Q_goal(1,1) = 0.01;
    Q_goal(2,2) = 0.01;
    
    l_obs = wobs*exp(- 0.5 * e_obs_pos' * wobs_width * e_obs_pos);
    l12obs_x = - wobs_width * l_obs * e_obs_pos;
    l12obs_xx = wobs_width^2 * l_obs*(e_obs_pos * e_obs_pos') - l_obs * wobs_width * (eye(2));
    
    l_obs_x = [l12obs_x' 0 0]';
    l_obs_xx = [l12obs_xx zeros(2) ; zeros(2,4)];
    
    % Ensure cost's positive definiteness
%     [V,D] = eig(l_obs_xx);
%     d = diag(D);
%     d(d<0) = 0;
%     l_obs_xx = V*diag(d)*V';
    
    l_sobs = 0;
    l_sobs_x = zeros(4,1);
    l_sobs_xx = zeros(4);
    
    for i = 1:n_sobs 
        e_sobs =  x_rob(1:2) - sobs_pos(:,i);
        Q_expected = (eye(2) - Q_sobs * sobs_Sigma(:,:,i))\Q_sobs;
        l_sobsi = wobs*exp(- 0.5 * e_sobs' *Q_expected * e_sobs);
        l12sobs_x = - l_sobsi * Q_expected * e_sobs;
        l12sobs_xx = l_sobsi * (Q_expected * e_sobs) * (Q_expected * e_sobs)' - l_sobsi * Q_expected;
        l_sobs = l_sobs + l_sobsi;
        l_sobs_x(1:2) = l_sobs_x(1:2) + l12sobs_x;
        l_sobs_xx(1:2,1:2) = l_sobs_xx(1:2,1:2) + l12sobs_xx;
    end
    
    % Ensure cost's positive definiteness
%     [V,D] = eig(l_sobs_xx);
%     d = diag(D);
%     d(d<0) = 0;
%     l_sobs_xx = V*diag(d)*V';


    l_goal = e_goal'*Q_goal*e_goal;
    l_goal_x = 2*Q_goal*e_goal;
    l_goal_xx = 2*Q_goal;

    
    l = l_goal + l_obs + l_sobs;
    l_x = [l_sobs_x; l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
    l_xx = [l_sobs_xx zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx + l_sobs_xx];
    
%     l_x = [zeros(4,1); l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
%     l_xx = [zeros(4) zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx + l_sobs_xx];
    
%     l_x = [zeros(4,1); l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
%     l_xx = [zeros(4) zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx l_sobs_xx ; zeros(4) l_sobs_xx l_obs_xx + l_sobs_xx];
% 
%     l = l_goal + l_obs + l_sobs;
%     l_x = [l_sobs_x; l_goal_x ; l_obs_x ];
%     l_xx = [l_sobs_xx zeros(4,8); zeros(4) l_goal_xx  zeros(4) ; zeros(4) zeros(4) l_obs_xx ];
    
    % Ensure cost's positive definiteness
    [V,D] = eig(l_xx);
    d = diag(D);
    d(d<0) = 0;
    l_xx = V*diag(d)*V';
else
    Q_goal(1,1) = 1;
    Q_goal(2,2) = 1;
    wobs = 0.1;
    l_obs = wobs*exp(- 0.5 * e_obs_pos' * wobs_width * e_obs_pos);
    l12obs_x = - wobs_width * l_obs * e_obs_pos;
    l12obs_xx = wobs_width^2 * l_obs*(e_obs_pos * e_obs_pos') - l_obs * wobs_width * (eye(2));
    
    l_obs_x = [l12obs_x' 0 0]';
    l_obs_xx = [l12obs_xx zeros(2) ; zeros(2,4)];
    
    % Ensure cost's positive definiteness
%     [V,D] = eig(l_obs_xx);
%     d = diag(D);
%     d(d<0) = 0;
%     l_obs_xx = V*diag(d)*V';
    
    l_sobs = 0;
    l_sobs_x = zeros(4,1);
    l_sobs_xx = zeros(4);
    
    for i = 1:n_sobs
        e_sobs =  x_rob(1:2) - sobs_pos(:,i);
        Q_expected = (eye(2) - Q_sobs * sobs_Sigma(:,:,i))\Q_sobs;
        l_sobsi = wobs*exp(- 0.5 * e_sobs' *Q_expected * e_sobs);
        l12sobs_x = - l_sobsi * Q_expected * e_sobs;
        l12sobs_xx = l_sobsi * (Q_expected * e_sobs) * (Q_expected * e_sobs)' - l_sobsi * Q_expected;
        l_sobs = l_sobs + l_sobsi;
        l_sobs_x(1:2) = l_sobs_x(1:2) + l12sobs_x;
        l_sobs_xx(1:2,1:2) = l_sobs_xx(1:2,1:2) + l12sobs_xx;
    end
    
    % Ensure cost's positive definiteness
%     [V,D] = eig(l_sobs_xx);
%     d = diag(D);
%     d(d<0) = 0;
%     l_sobs_xx = V*diag(d)*V';


    l_goal = e_goal'*Q_goal*e_goal;
    l_goal_x = 2*Q_goal*e_goal;
    l_goal_xx = 2*Q_goal;

    
    l = l_goal + l_obs + l_sobs;
    l_x = [l_sobs_x; l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
    l_xx = [l_sobs_xx zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx + l_sobs_xx];
    
    
%     l_x = [zeros(4,1); l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
%     l_xx = [zeros(4) zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx zeros(4) ; zeros(4) zeros(4) l_obs_xx + l_sobs_xx];
    
%     l_x = [zeros(4,1); l_goal_x - l_sobs_x; l_obs_x - l_sobs_x];
%     l_xx = [zeros(4) zeros(4,8); zeros(4) l_goal_xx + l_sobs_xx l_sobs_xx ; zeros(4) l_sobs_xx l_obs_xx + l_sobs_xx];
%     
%     l = l_goal + l_obs + l_sobs;
%     l_x = [l_sobs_x; l_goal_x ; l_obs_x ];
%     l_xx = [l_sobs_xx zeros(4,8); zeros(4) l_goal_xx  zeros(4) ; zeros(4) zeros(4) l_obs_xx ];
    
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

