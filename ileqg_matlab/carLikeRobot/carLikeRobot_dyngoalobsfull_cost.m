function [ l, l_x, l_xx, l_u, l_uu, l_ux  ] = carLikeRobot_dyngoalobsfull_cost( x, u, t, sobs_pos, sobs_Sigma)
% CARLIKEROBOT_COST 
% Cost penalizing quadratic measure 

sx = length(x);
su = length(u);

x_rob = x(1:5); 
e_goal = x(6:10);
e_obs_pos = x(11:12);

Q_goal = zeros(5);

Q_goal(1,1) = 1e1;
Q_goal(2,2) = 1e1;
wobs = 0.1;
wobs_width = 1/0.005;

% Static obstacles
Q_sobs = 1*eye(2);

if(isnan(t))
    Q_sobs = Q_sobs*0.01;
    Q_goal = Q_goal*0.01;
    wobs = wobs*0.01;
end

% Exponential functional to dynamic obstacle obs
l_obs = wobs*exp(- 0.5 * e_obs_pos' * wobs_width * e_obs_pos);
l12obs_x = - wobs_width * l_obs * e_obs_pos;
l12obs_xx = wobs_width * l_obs*(e_obs_pos * e_obs_pos') - l_obs * wobs_width * eye*(2);
l_obs_x = [l12obs_x' 0 0]';
l_obs_xx = [l12obs_xx zeros(2) ; zeros(2,4)];

% Static obstacles sobs
n_sobs = size(sobs_pos, 2);
l_sobs = 0;
l_sobs_x = zeros(5,1);
l_sobs_xx = zeros(5);

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

% Quadratic goal
l_goal = e_goal'*Q_goal*e_goal;
l_goal_x = 2*Q_goal*e_goal;
l_goal_xx = 2*Q_goal;


% Build up cost
l = l_goal + l_obs + l_sobs;
l_obs_x(1:2) = l_obs_x(1:2) - l_sobs_x(1:2);
l_x = [zeros(5,1); l_goal_x - l_sobs_x; l_obs_x];
l_obs_xx(1:2,1:2) = l_obs_xx(1:2,1:2) + l_sobs_xx(1:2,1:2);
l_xx = [zeros(5) zeros(5,9); zeros(5) l_goal_xx + l_sobs_xx zeros(5,4) ; zeros(4,5) zeros(4,5) l_obs_xx];

% Ensure cost's positive definiteness
[V,D] = eig(l_xx);
d = diag(D);
d(d<1e-6) = 1e-6;
l_xx = V*diag(d)*V';


% Control costs
if ~isnan(t)
    R = 1e-3 * eye(su);
    
    l_control = u'*R*u;
    
    l = l + l_control;

    l_u = 2*R*u;

    l_uu = 2*R;

    l_ux = zeros(su, sx);
end


end
