close all;
startup;

%% Optimization parameters
T = 200;
dt = 0.01;
u0 = zeros(2,1);
max_iter = 1000;
uMin = -Inf;
uMax = Inf;

% Initial states
x_rob = [0 0 0 0]';
x_goal = [0.1 0.1 1 -0.5]';
x_obs = [0.15 0.15 1 -0.5]';
x0 = [x_rob ; (x_goal - x_rob) ; (x_obs - x_rob)];

% Static obstacles: position and Sigma
sobs_pos = [0.8 -0.3 ; 0.8 -0.6]';
sobs_Sigma(:,:,1) =  0.001*eye(2);
sobs_Sigma(:,:,2) =  0.005*eye(2);

% Cost function
costWithObstacles = ...
    @(x, u , t) pointMass_dyngoalobsfull_cost(x, u, t, sobs_pos, sobs_Sigma);
costWithoutObstacles = ... % For initialization
    @(x, u , t) pointMass_dyngoalobsfull_cost(x, u, t, [], []); 
                                      
% 4 different dynamics with different noise levels
Gamma1 = zeros(length(x0));
Gamma2 = zeros(length(x0));
Gamma1(9:12,9:12) = 0.8*eye(4);
Gamma2(5:8,5:8) = 0.8*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;
dynamicsGamma1 = @(x, u) pointMass_dyngoalobsfullvar_dyn(x, u, Gamma);
Gamma1(9:12,9:12) = 0.5*eye(4);
Gamma2(5:8,5:8) = 0.5*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;
dynamicsGamma2 = @(x, u) pointMass_dyngoalobsfullvar_dyn(x, u, Gamma);
Gamma1(9:12,9:12) = 0.1*eye(4);
Gamma2(5:8,5:8) = 0.2*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;
dynamicsGamma3 = @(x, u) pointMass_dyngoalobsfullvar_dyn(x, u, Gamma);
Gamma1(9:12,9:12) = 0.01*eye(4);
Gamma2(5:8,5:8) = 0.01*eye(4);
Gamma(:,:,1) = Gamma1;
Gamma(:,:,2) = Gamma2;
dynamicsGamma4 = @(x, u) pointMass_dyngoalobsfullvar_dyn(x, u, Gamma);

% Plot functions
plotWithoutObstacles = ...
    @(x, u , L, varargin) pointMass_dyngoalobsfull_plot(x, u, L, ...
                                                        [], [], varargin{:});
plotWithObstacles = ...
    @(x, u , L, varargin) pointMass_dyngoalobsfull_plot(x, u, L, ...
                                          sobs_pos, sobs_Sigma, varargin{:});

% Define LQ solver
gamma_neg = 2;
gamma_pos = 4;
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
kccOptMeanNVarPVar = @(lqProb, varargin) ...
          kcc_multi([1 gamma_pos^2 ; 1 -(gamma_neg^2)], lqProb, varargin{:});

%% Solve
% Initialize trajectory with expected solver in scene without obstacles
[x_init, u_init, L_init, cost_init, eps_init] = ...
            iterativeLQSolver(dynamicsGamma1, costWithoutObstacles, ...
                              ExpectedCost, dt, T, x0, u0, uMin, uMax, ...
                              max_iter, plotWithoutObstacles);

% Solve with noise level 1
[x_mvnp1, u_mvnp1, L_mvnp1, cost_mvnp1, eps_mvnp1] = ...
            iterativeLQSolver(dynamicsGamma1, costWithObstacles, ...
                              kccOptMeanNVarPVar, dt, T, x0, u_init, ...
                              uMin, uMax, max_iter, plotWithObstacles);

% Solve with noise level 2
[x_mvnp2, u_mvnp2, L_mvnp2, cost_mvnp2, eps_mvnp2] = ...
            iterativeLQSolver(dynamicsGamma2, costWithObstacles, ...
                              kccOptMeanNVarPVar, dt, T, x0, u_init, ...
                              uMin, uMax, max_iter, plotWithObstacles);

% Solve with noise level 3
[x_mvnp3, u_mvnp3, L_mvnp3, cost_mvnp3, eps_mvnp3] = ...
            iterativeLQSolver(dynamicsGamma3, costWithObstacles, ...
                              kccOptMeanNVarPVar, dt, T, x0, u_init, ...
                              uMin, uMax, max_iter, plotWithObstacles);
                          
% Solve with noise level 4
[x_mvnp4, u_mvnp4, L_mvnp4, cost_mvnp4, eps_mvnp4] = ...
            iterativeLQSolver(dynamicsGamma4, costWithObstacles, ...
                              kccOptMeanNVarPVar, dt, T, x0, u_init, ...
                              uMin, uMax, max_iter, plotWithObstacles);

%% Simulate solutions and get statistics
n_sims = 500;
idx_rob = 1;
idx_des = 5;
idx_obs = 9;

% Trajectories and cost stats
[x_trajs_mvnp1, u_trajs_mvnp1, sobs_mvnp1] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvnp1, u_mvnp1, L_mvnp1, ...
                             dynamicsGamma1, sobs_pos, sobs_Sigma);
[cost_stats_mvnp1] = ...
    compute_cost_stats(costWithObstacles, x_trajs_mvnp1, u_trajs_mvnp1, dt);


[x_trajs_mvnp2, u_trajs_mvnp2, sobs_mvnp2] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvnp2, u_mvnp2, L_mvnp2, ...
                             dynamicsGamma2, sobs_pos, sobs_Sigma);
[cost_stats_mvnp2] = ...
    compute_cost_stats(costWithObstacles, x_trajs_mvnp2, u_trajs_mvnp2, dt);


[x_trajs_mvnp3, u_trajs_mvnp3, sobs_mvnp3] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvnp3, u_mvnp3, L_mvnp3, ...
                             dynamicsGamma3, sobs_pos, sobs_Sigma);
[cost_stats_mvnp3] = ...
    compute_cost_stats(costWithObstacles, x_trajs_mvnp3, u_trajs_mvnp3, dt);


[x_trajs_mvnp4, u_trajs_mvnp4, sobs_mvnp4] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvnp4, u_mvnp4, L_mvnp4, ...
                             dynamicsGamma4, sobs_pos, sobs_Sigma);
[cost_stats_mvnp4] = ...
    compute_cost_stats(costWithObstacles, x_trajs_mvnp4, u_trajs_mvnp4, dt);

% Collision stats
rob_size = 0.05;
[avg_colls_mvnp1, var_colls_mvnp1] = ...
    compute_collisions(x_trajs_mvnp1, idx_obs, rob_size);
[avg_colls_mvnp2, var_colls_mvnp2] = ...
    compute_collisions(x_trajs_mvnp2, idx_obs, rob_size);
[avg_colls_mvnp3, var_colls_mvnp3] = ...
    compute_collisions(x_trajs_mvnp3, idx_obs, rob_size);
[avg_colls_mvnp4, var_colls_mvnp4] = ...
    compute_collisions(x_trajs_mvnp4, idx_obs, rob_size);


%% PLOT 4 solutions
xmin = -0.1;
xmax = 1.5;
ymin = -0.85;
ymax = 0.35;

transparency = 0.5;
color_rob = [[178,223,138]./255 transparency];
color_obs = [[251,154,153]./255 transparency];
color_goal = [[166,206,227]./255 transparency];
color_sobs =  [[253,191,111]./255 transparency];


fig_pointmass_full  = figure('units','normalized', ...
                             'outerposition',[0 0 1 0.4]);

ax1 = axes('Parent',fig_pointmass_full, ...
           'Position',[0.0502059202059202 0.239401496259352 0.2 0.7],...
           'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
           'FontSize',20,...
           'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_mvnp1, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvnp1, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_mvnp1, color_sobs);
plot_n_simulations_idx(x_trajs_mvnp1, idx_rob, [], color_rob);
plotWithObstacles(x_mvnp1, u_mvnp1, L_mvnp1, ax1);

goal = findobj('Tag','goal');
obs = findobj('Tag','dynamic obstacle');
robot = findobj('Tag','robot');
gain_obs = findobj('Tag','gain obstacle');
gain_goal = findobj('Tag','gain goal');
static_obs = findobj('Tag','static obstacle');

xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
set(ax1,'YTick',[-0.6,-0.3,0,0.3]);
set(ax1,'XTick',[0,0.3,0.6,0.9,1.2,1.5]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
ylabel(ax1, 'ay1');

ax2 = axes('Parent',fig_pointmass_full, ...
    'Position',[0.279488964648537 0.239401496259352 0.2 0.7],...
    'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
    'FontSize',20,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_mvnp2, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvnp2, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_mvnp2, color_sobs);
plot_n_simulations_idx(x_trajs_mvnp2, idx_rob, [], color_rob);
plotWithObstacles(x_mvnp2, u_mvnp2, L_mvnp2, ax2);

xlim(ax2,[xmin xmax]);
ylim(ax2,[ymin ymax]);
set(ax2,'YTick',[-0.6,-0.3,0,0.3]);
set(ax2,'XTick',[0,0.3,0.6,0.9,1.2,1.5]);
box(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'ax2');

ax3 = axes('Parent',fig_pointmass_full,...
    'Position',[0.509652509652508 0.239401496259352 0.2 0.7],...
    'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
    'FontSize',20,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_mvnp3, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvnp3, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_mvnp3, color_sobs);
plot_n_simulations_idx(x_trajs_mvnp3, idx_rob, [], color_rob);
plotWithObstacles(x_mvnp3, u_mvnp3, L_mvnp3, ax3);

xlim(ax3,[xmin xmax]);
ylim(ax3,[ymin ymax]);
set(ax3,'YTick',[-0.6,-0.3,0,0.3]);
set(ax3,'XTick',[0,0.3,0.6,0.9,1.2,1.5]);
box(ax3, 'on');
grid(ax3, 'on');
xlabel(ax3, 'ax3');

ax4 = axes('Parent',fig_pointmass_full,...
    'Position',[0.74002574002574 0.239401496259352 0.2 0.7],...
    'PlotBoxAspectRatio',[1.33333333333333 1 1.66666666666667],...
    'FontSize',20,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_mvnp4, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvnp4, idx_obs, idx_rob, color_obs);
plot_n_static_obs(sobs_mvnp4, color_sobs);
plot_n_simulations_idx(x_trajs_mvnp4, idx_rob, [], color_rob);
plotWithObstacles(x_mvnp4, u_mvnp4, L_mvnp4, ax4);

xlim(ax4,[xmin xmax]);
ylim(ax4,[ymin ymax]);
set(ax4,'YTick',[-0.6,-0.3,0,0.3]);
set(ax4,'XTick',[0,0.3,0.6,0.9,1.2,1.5]);
box(ax4, 'on');
grid(ax4, 'on');
xlabel(ax4, 'ax4');

legend1 = legend([robot(1) goal(1) obs(1) static_obs(1)], ...
  {'robot','goal', 'dynamicobstacle', 'staticobstaclessssss'});
set(legend1,'Location','SouthWest','FontSize',12);

file = '../../repositories/dhri/drafts/onILEQR/figures/simulation_pointMassFull_varnoise.eps';
set(fig_pointmass_full,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_pointmass_full, '-depsc','-painters', file);