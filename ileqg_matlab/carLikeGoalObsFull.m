close all;

%% Optimization parameters
T = 200;
dt = 0.01;
theta = 1;
u0 = zeros(2,1);
max_iter = 1000;
uMin = -Inf;
uMax = Inf;

% Initial state
x_rob = [0 0 0 0.01 0.01]';
x_rob_taskspace = [x_rob(1) x_rob(2) ...
                   x_rob(4)*cos(x_rob(3)) x_rob(4)*sin(x_rob(3))]';
x_goal = [0.01 0.01 0.1 0.4 0.8]';
x_obs = [0.05 0.15 0.25 0]';
x0 = [x_rob ; (x_goal - x_rob); (x_obs - x_rob_taskspace)];

% Static obstacles: position and uncertainty
sobs_pos = [0.6 0.3]';
sobs_Sigma(:,:,1) =  0.0001*eye(2);

% Cost function
costWithObstacles = ...
    @(x, u , t) carLikeRobot_dyngoalobsfull_cost(x, u, t, ...
                                                 sobs_pos, sobs_Sigma);
costWithoutObstacles = ... % for initialization
    @(x, u , t) carLikeRobot_dyngoalobsfull_cost(x, u, t, [], []);

% Dynamic function
dynamics_function = @carLikeRobot_dyngoalobsfull_dyn;

% Plot functions
plotWithoutObstacles = ...
    @(x, u , L, varargin) carLikeRobot_dyngoalobsfull_plot(x, u, L, ...
                                                          [], [], varargin{:});
plotWithObstacles = ...
    @(x, u , L, varargin) carLikeRobot_dyngoalobsfull_plot(x, u, L, ...
                                            sobs_pos, sobs_Sigma, varargin{:});

% Define LQ solvers
gamma_neg = 3;
gamma_pos = 200;
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
kccOptMeanNVarPVar = ...
    @(lqProb, varargin) kcc_multi([1 gamma_pos^2 ; 1 -(gamma_neg^2)], ...
                                  lqProb, varargin{:});
riskSeekingAverse = ...
    @(lqProb, varargin)leqr_multi([0.001 ; -1000], lqProb, varargin{:});


%% Solve
% Initialization
[x_exp, u_exp, L_exp, cost_exp, eps_exp] = ...
    iterativeLQSolver(dynamics_function, costWithObstacles, ExpectedCost, ...
                      dt, T, x0, u0, uMin, uMax, max_iter, plotWithObstacles);        
% Solution
[x_mvnp, u_mvnp, L_mvnp, cost_mvnp, eps_mvnp] = ...
    iterativeLQSolver(dynamics_function, costWithObstacles, ...
                      kccOptMeanNVarPVar, dt, T, x0, u_exp, uMin, uMax, ...
                      max_iter, plotWithObstacles);


%% Simulations
n_sims = 100;
idx_rob = 1;
idx_des = 6;
idx_obs = 11;
transparency = 0.5;
color_rob = [[178,223,138]./255 transparency];
color_obs = [[251,154,153]./255 transparency];
color_goal = [[166,206,227]./255 transparency];
color_sobs =  [[253,191,111]./255 transparency];

[x_trajs_mvnp, u_trajs_mvnp, sobs_mvnp] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvnp, u_mvnp, L_mvnp, ...
                             dynamics_function, sobs_pos, sobs_Sigma);
[cost_stats_mvnp] = ...
    compute_cost_stats(costWithObstacles, x_trajs_mvnp, u_trajs_mvnp, dt);

[x_trajs_exp, u_trajs_exp, sobs_exp] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_exp, u_exp, L_exp, ...
                             dynamics_function, sobs_pos, sobs_Sigma);
[cost_stats_exp] = ...
    compute_cost_stats(costWithObstacles, x_trajs_exp, u_trajs_exp, dt);

% Collision checking
rob_size = 0.05;
[avg_colls_mvnp, var_colls_mvnp] = ...
    compute_collisions(x_trajs_mvnp, idx_obs, rob_size);
[avg_colls_exp, var_colls_exp] = ...
    compute_collisions(x_trajs_exp, idx_obs, rob_size);

%% PLOT
xmin = -0.1;
xmax = 0.9;
ymin = -0.1;
ymax = 0.5;

fig_pointmass_full  = figure('units','normalized', ...
                             'outerposition',[0 0 0.25 0.65]);

ax1 = axes('Parent',fig_pointmass_full, ...
    'Position',[0.13573407 0.55718206 0.8407202216 0.43950701],...
    'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
    'FontSize',22,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_exp, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_exp, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_exp, color_sobs);
plot_n_simulations_idx(x_trajs_exp, idx_rob, [], color_rob);
plotWithObstacles(x_exp, u_exp, L_exp, ax1);
xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
ylabel(ax1, 'ay1');
set(ax1,'YTick',[0,0.2,0.4,0.6]);
set(ax1,'XTick',[0,0.2,0.4,0.6,0.8]);
title('t1');


ax2 = axes('Parent',fig_pointmass_full, ...
    'Position',[0.13573407202 0.011235132662 0.84072022160 0.44530443808],...
    'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
    'FontSize',22,...
    'DataAspectRatio',[1 1 1]);

plot_n_simulations_idx(x_trajs_mvnp, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvnp, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_mvnp, color_sobs);
plot_n_simulations_idx(x_trajs_mvnp, idx_rob, [], color_rob);
plotWithObstacles(x_mvnp, u_mvnp, L_mvnp, ax2);
xlim(ax2,[xmin xmax]);
ylim(ax2,[ymin ymax]);
box(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'ax2');
ylabel(ax2, 'ay2');
set(ax2,'YTick',[0,0.2,0.4,0.6]);
set(ax2,'XTick',[0,0.2,0.4,0.6,0.8]);
title('t2');

goal = findobj('Tag','goal');
obs = findobj('Tag','dynamic obstacle');
robot = findobj('Tag','robot');
static_obs = findobj('Tag','static obstacle');

legend1 = legend([robot(1) goal(1) obs(1) static_obs(1)], ...
                 {'robot','goal', 'dynamicobstacle',  'staticobstaclesssss'});
set(legend1,'Location','NorthWest','FontSize',14);

file = '../../../repositories/dhri/drafts/onILEQR/figures/simulation_carLikeFull.eps';
set(fig_pointmass_full,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_pointmass_full, '-depsc','-painters', file);