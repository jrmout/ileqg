close all;
startup;

%% Optimization parameters
T = 50;
dt = 0.01;
u0 = zeros(2,1);
max_iter = 1000;
uMin = -Inf;
uMax = Inf;

x_rob = [0 0 0 0]';
x_goal = [0.01 0.01 1 -0.5]';
x_obs = [0.02 0.02 1 -0.5]';
x0 = [x_rob ; (x_goal - x_rob) ; (x_obs - x_rob)];

gamma_neg = 5;
gamma_pos = 5;

% Cost function
cost_function = @pointMass_dyngoalobs_cost;
% Dynamics function
dynamics_function = @pointMass_dyngoalobs_dyn_2noises;

% Define LQ solvers
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
kccOptMeanNVarPVar = ...
    @(lqProb, varargin)kcc_multi([1 2*gamma_pos^2 ; 1 -2*(gamma_neg^2)], ...
                                 lqProb, varargin{:});
riskSeekingAverse = ...
    @(lqProb, varargin)leqr_multi([3e2 ; -1e2], lqProb, varargin{:});

%% Solve
[x_exp, u_exp, L_exp, cost_exp, eps_exp] = ...
    iterativeLQSolver(dynamics_function, cost_function, ExpectedCost, ...
                      dt, T, x0, u0, uMin, uMax, max_iter);
[x_mv, u_mv, L_mv, cost_mv, eps_mv] = ...
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNVarPVar, ...
                      dt, T, x0, u_exp, uMin, uMax, max_iter);
[x_r, u_r, L_r, cost_r, eps_r] = ...
    iterativeLQSolver(dynamics_function, cost_function, riskSeekingAverse, ...
                      dt, T, x0, u_exp, uMin, uMax, max_iter);

%% Simulations
n_sims = 500;
idx_rob = 1;
idx_des = 5;
idx_obs = 9;

[x_trajs_mv, u_trajs_mv] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mv, u_mv, ...
                             L_mv, dynamics_function);
[cost_stats_mv] = ...
    compute_cost_stats(cost_function, x_trajs_mv, u_trajs_mv, dt);

[x_trajs_r, u_trajs_r] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_r, u_r, ...
                             L_r, dynamics_function);
[cost_stats_r] = ...
    compute_cost_stats(cost_function, x_trajs_r, u_trajs_r, dt);

[x_trajs_exp, u_trajs_exp] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_exp, u_exp, ...
                             L_exp, dynamics_function);
[cost_stats_exp] = ...
    compute_cost_stats(cost_function, x_trajs_exp, u_trajs_exp, dt);

% Collision checking
rob_size = 0.02;
[avg_colls_mv, var_colls_mv] = ...
    compute_collisions(x_trajs_mv, idx_obs, rob_size);
[avg_colls_r, var_colls_r] = ...
    compute_collisions(x_trajs_r, idx_obs, rob_size);
[avg_colls_exp, var_colls_exp] = ...
    compute_collisions(x_trajs_exp, idx_obs, rob_size);

%% PLOT
xmin = -0.05;
xmax = 0.5;
ymin = -0.275;
ymax = 0.05;

transparency = 0.5;
color_rob = [[178,223,138]./255 transparency];
color_obs = [[251,154,153]./255 transparency];
color_goal = [[166,206,227]./255 transparency];

fig_dyngoal_2s  = figure('units','normalized','outerposition',[0 0 0.3 0.45]);

ax1 = axes('Parent',fig_dyngoal_2s, ...
    'Position',[0.06577405121 0.5451615651 0.4382443 0.4681139755],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'FontSize',16,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_exp, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_exp, idx_obs, idx_rob, color_obs);
plot_n_simulations_idx(x_trajs_exp, idx_rob, [], color_rob);
pointMass_dyngoalobs_plot(x_exp, u_exp, L_exp, ax1);

goal = findobj('Tag','goal');
obs = findobj('Tag','dynamic obstacle');
robot = findobj('Tag','robot');
gain_obs = findobj('Tag','gain obstacle');
gain_goal = findobj('Tag','gain goal');

xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
ylabel(ax1, 'ay1');
legend1 = legend([robot(1) goal(1) obs(1) ], ...
                 {'robot','goal', 'dynamicobstacle'});
legend('boxoff');
set(legend1,'Position',[0.12 0.139650 0.13839 0.316084],'FontSize',16);

ax2 = axes('Parent',fig_dyngoal_2s,'YTickLabel',{'','','','',''}, ...
    'Position',[0.535017221584387 0.544226521575613 0.4382443184773 0.468113975576625],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'FontSize',16,...
    'DataAspectRatio',[1 1 1]);

plot_n_simulations_idx(x_trajs_mv, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mv, idx_obs, idx_rob, color_obs);
plot_n_simulations_idx(x_trajs_mv, idx_rob, [], color_rob);
pointMass_dyngoalobs_plot(x_mv, u_mv, L_mv, ax2);
xlim(ax2,[xmin xmax]);
ylim(ax2,[ymin ymax]);
set(ax2,'YTickLabel',{'','','','',''})
box(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'ax2');

ax3 = axes('Parent',fig_dyngoal_2s, ...
    'Position',[0.5354707233 0.05327670432 0.43824431 0.4681139755],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'FontSize',16,...
    'DataAspectRatio',[1 1 1]);
plot_n_simulations_idx(x_trajs_r, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_r, idx_obs, idx_rob, color_obs);
plot_n_simulations_idx(x_trajs_r, idx_rob, [], color_rob);
pointMass_dyngoalobs_plot(x_r, u_r, L_r, ax3);
xlim(ax3,[xmin xmax]);
ylim(ax3,[ymin ymax]);
box(ax3, 'on');
grid(ax3, 'on');
xlabel(ax3, 'ax3');
ylabel(ax3, 'ay3');

file = '../../repositories/dhri/drafts/onILEQR/figures/simulation_dyngoalobs_2s.eps';
set(fig_dyngoal_2s,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_dyngoal_2s, '-depsc','-painters', file);