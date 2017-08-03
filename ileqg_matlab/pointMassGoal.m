close all;
startup

% Optimization parameters
T = 50;
dt = 0.01;
theta = 1;
x0 = [0 0 0 0 0.1 0.1 1 -0.5]';
u0 = zeros(2,1);
max_iter = 1000;
uMin = -Inf;
uMax = Inf;

gamma = 9;

% Dynamics
dynamics_function = @pointMass_dyngoal_dyn;

% Cost
cost_function = @pointMass_dyngoal_cost;

% Plot (during optimization)
plot_function = @pointMass_dyngoal_plot;

% Define LQ solvers
% E[J]
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
% E[J]-(gamma_neg)/2 Var[J]
kccOptMeanNVar = @(lqProb, varargin)kcc_multi(... 
                                [1 -(gamma)/2], lqProb, varargin{:});
% E[J]+(gamma_pos)/2 Var[J]
kccOptMeanPVar = @(lqProb, varargin)kcc_multi(... 
                                [1 gamma /2], lqProb, varargin{:});
% E[J]-(gamma_neg^2)/6 k3[J]
kccOptMeanNThird = @(lqProb, varargin)kcc_multi(... 
                                [1 0 -(gamma^2) / 6], lqProb, varargin{:});
% E[J]+(gamma_pos^2)/6 k3[J]
kccOptMeanPThird = @(lqProb, varargin)kcc_multi(... 
                                [1 0 gamma^2 / 6], lqProb, varargin{:});
% E[J]-(gamma_neg^3)/24 k4[J]
kccOptMeanNFourth = @(lqProb, varargin)kcc_multi(...
                                [1 0 0 -(gamma^3) / 24], lqProb, varargin{:});
% E[J]+(gamma_pos^3)/24 k4[J]
kccOptMeanPFourth = @(lqProb, varargin)kcc_multi(...
                                [1 0 0 gamma^3 / 24], lqProb, varargin{:});
% E[exp(theta_pos/2 J)]
riskAverse = @(lqProb, varargin)leqr_multi([gamma], lqProb, varargin{:});
% E[exp(theta_neg/2 J)]
riskSeeking = @(lqProb, varargin)leqr_multi([-gamma], lqProb, varargin{:});

%% Solve iteratively
[x_exp, u_exp, L_exp, cost_exp, eps_exp] = ... 
    iterativeLQSolver(dynamics_function, cost_function, ExpectedCost, ...
                   dt, T, x0, u0, uMin, uMax, max_iter, plot_function);
[x_mvn, u_mvn, L_mvn, cost_mvn, eps_mvn] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNVar, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_mvp, u_mvp, L_mvp, cost_mvp, eps_mvp] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanPVar, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_m3n, u_m3n, L_m3n, cost_m3n, eps_m3n] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNThird, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_m3p, u_m3p, L_m3p, cost_m3p, eps_m3p] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanPThird, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_m4n, u_m4n, L_m4n, cost_m4n, eps_m4n] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNFourth, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_m4p, u_m4p, L_m4p, cost_m4p, eps_m4p] = ... 
    iterativeLQSolver(dynamics_function, cost_function, kccOptMeanPFourth, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_ra, u_ra, L_ra, cost_ra, eps_ra] = ... 
    iterativeLQSolver(dynamics_function, cost_function, riskAverse, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
[x_rs, u_rs, L_rs, cost_rs, eps_rs] = ... 
    iterativeLQSolver(dynamics_function, cost_function, riskSeeking, ...
                   dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);

%% Simulations
n_sims = 500;
idx_rob = 1;
idx_des = 5;

[x_trajs_exp, u_trajs_exp] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_exp, u_exp, ...
                             L_exp, dynamics_function);
[cost_stats_exp] = ...
    compute_cost_stats(cost_function, x_trajs_exp, u_trajs_exp, dt);

[x_trajs_mvn, u_trajs_mvn] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvn, u_mvn, ...
                             L_mvn, dynamics_function);
[cost_stats_mvn] = ...
    compute_cost_stats(cost_function, x_trajs_mvn, u_trajs_mvn, dt);

[x_trajs_mvp, u_trajs_mvp] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_mvp, u_mvp, ...
                             L_mvp, dynamics_function);
[cost_stats_mvp] = ...
    compute_cost_stats(cost_function, x_trajs_mvp, u_trajs_mvp, dt);

[x_trajs_m3n, u_trajs_m3n] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_m3n, u_m3n, ...
                             L_m3n, dynamics_function);
[cost_stats_m3n] = ...
    compute_cost_stats(cost_function, x_trajs_m3n, u_trajs_m3n, dt);

[x_trajs_m3p, u_trajs_m3p] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_m3p, u_m3p, ...
                             L_m3p, dynamics_function);
[cost_stats_m3p] = ...
    compute_cost_stats(cost_function, x_trajs_m3p, u_trajs_m3p, dt);

[x_trajs_m4n, u_trajs_m4n] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_m4n, u_m4n, ...
                             L_m4n, dynamics_function);
[cost_stats_m4n] = ...
    compute_cost_stats(cost_function, x_trajs_m4n, u_trajs_m4n, dt);

[x_trajs_m4p, u_trajs_m4p] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_m4p, u_m4p, ...
                             L_m4p, dynamics_function);
[cost_stats_m4p] = ...
    compute_cost_stats(cost_function, x_trajs_m4p, u_trajs_m4p, dt);

[x_trajs_ra, u_trajs_ra] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_ra, u_ra, ...
                             L_ra, dynamics_function);
[cost_stats_ra] = ...
    compute_cost_stats(cost_function, x_trajs_ra, u_trajs_ra, dt);

[x_trajs_rs, u_trajs_rs] = ...
    stochastic_n_simulations(n_sims, dt, T, x0, x_rs, u_rs, ...
                             L_rs, dynamics_function);
[cost_stats_rs] = ...
    compute_cost_stats(cost_function, x_trajs_rs, u_trajs_rs, dt);

figure_cost_stas = figure;
subplot(2,2,1);
bar([cost_stats_exp(1),cost_stats_mvn(1),cost_stats_mvp(1), ...
     cost_stats_m3n(1),cost_stats_m3p(1),cost_stats_m4n(1), ...
     cost_stats_m4p(1),cost_stats_rs(1),cost_stats_ra(1)]);
title('E[J]');
subplot(2,2,2);
bar([cost_stats_exp(2),cost_stats_mvn(2),cost_stats_mvp(2), ...
     cost_stats_m3n(2),cost_stats_m3p(2),cost_stats_m4n(2), ...
     cost_stats_m4p(2),cost_stats_rs(2),cost_stats_ra(2)]);
title('Var[J]');
subplot(2,2,3);
bar([cost_stats_exp(3),cost_stats_mvn(3),cost_stats_mvp(3), ...
     cost_stats_m3n(3),cost_stats_m3p(3),cost_stats_m4n(3), ...
     cost_stats_m4p(3),cost_stats_rs(3),cost_stats_ra(3)]);
title('k3[J]');
subplot(2,2,4);
bar([cost_stats_exp(4),cost_stats_mvn(4),cost_stats_mvp(4), ...
     cost_stats_m3n(4),cost_stats_m3p(4),cost_stats_m4n(4), ...
     cost_stats_m4p(4),cost_stats_rs(4),cost_stats_ra(4)]);
title('k4[J]');
 
%% PLOT
xmin = -0.1;
xmax = 0.65;
ymin = -0.3;
ymax = 0.15;

transparency = 0.5;
color_rob = [[178,223,138]./255 transparency];
color_goal = [[166,206,227]./255 transparency];

fig_dyngoal  = figure('units','normalized','outerposition',[0 0 0.5 0.3]);

ax1 = axes('Parent',fig_dyngoal,'Position',[0.057 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_exp, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_exp, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_exp, u_exp, L_exp, ax1);
xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
goal = findobj('Tag','goal');
robot = findobj('Tag','robot');
legend1 = legend(ax1,'show', [robot(1) goal(1)], {'l11111111111','l2222222222'});
legend('boxoff');
set(legend1,'Position',[0.04043 0.13965 0.13839 0.31608],'FontSize',13);


ax2 = axes('Parent',fig_dyngoal,'Position',[0.245 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_mvn, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvn, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_mvn, u_mvn, L_mvn, ax2);
xlim(ax2,[xmin xmax]);
ylim(ax2,[ymin ymax]);
set(ax2,'YTickLabel',{'','','','',''})
box(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'ax2');

ax3 = axes('Parent',fig_dyngoal,'Position',[0.245 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_mvp, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_mvp, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_mvp, u_mvp, L_mvp, ax3);
xlim(ax3,[xmin xmax]);
ylim(ax3,[ymin ymax]);
box(ax3, 'on');
grid(ax3, 'on');
xlabel(ax3, 'ax3');

ax4 = axes('Parent',fig_dyngoal,'Position',[0.439 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_m3n, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_m3n, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_m3n, u_m3n, L_m3n, ax4);
set(ax4,'YTickLabel',{'','','','',''})
xlim(ax4,[xmin xmax]);
ylim(ax4,[ymin ymax]);
box(ax4, 'on');
grid(ax4, 'on');
xlabel(ax4, 'ax4');

ax5 = axes('Parent',fig_dyngoal,'Position',[0.439 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_m3p, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_m3p, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_m3p, u_m3p, L_m3p, ax5);
set(ax5,'YTickLabel',{'','','','',''})
xlim(ax5,[xmin xmax]);
ylim(ax5,[ymin ymax]);
box(ax5, 'on');
grid(ax5, 'on');
xlabel(ax5, 'ax5');

ax6 = axes('Parent',fig_dyngoal,'Position',[0.63 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_m4n, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_m4n, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_m4n, u_m4n, L_m4n, ax6);
set(ax6,'YTickLabel',{'','','','',''})
xlim(ax6,[xmin xmax]);
ylim(ax6,[ymin ymax]);
box(ax6, 'on');
grid(ax6, 'on');
xlabel(ax6, 'ax6');

ax7 = axes('Parent',fig_dyngoal,'Position',[0.63 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_m4p, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_m4p, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_m4p, u_m4p, L_m4p, ax7);
set(ax7,'YTickLabel',{'','','','',''})
xlim(ax7,[xmin xmax]);
ylim(ax7,[ymin ymax]);
box(ax7, 'on');
grid(ax7, 'on');
xlabel(ax7, 'ax7');

ax8 = axes('Parent',fig_dyngoal,'Position',[0.82 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_rs, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_rs, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_rs, u_rs, L_rs, ax8);
set(ax8,'YTickLabel',{'','','','',''})
xlim(ax8,[xmin xmax]);
ylim(ax8,[ymin ymax]);
box(ax8, 'on');
grid(ax8, 'on');
xlabel(ax8, 'ax8');

ax9 = axes('Parent',fig_dyngoal,'Position',[0.82 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
plot_n_simulations_idx(x_trajs_ra, idx_des, idx_rob, color_goal);
plot_n_simulations_idx(x_trajs_ra, idx_rob, [], color_rob);
pointMass_dyngoal_plot(x_ra, u_ra, L_ra, ax9);
set(ax9,'YTickLabel',{'','','','',''})
xlim(ax9,[xmin xmax]);
ylim(ax9,[ymin ymax]);
box(ax9, 'on');
grid(ax9, 'on');
xlabel(ax9, 'ax9');


file = '../../dhri/drafts/onILEQR/figures/simulation_dyngoal.eps';
set(fig_dyngoal,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_dyngoal, '-depsc2','-painters', file);


figure_cost_stas = figure;
subplot(2,2,1);
bar([cost_stats_exp(1),cost_stats_mvn(1),cost_stats_mvp(1), ...
     cost_stats_m3n(1),cost_stats_m3p(1),cost_stats_m4n(1), ...
     cost_stats_m4p(1),cost_stats_rs(1),cost_stats_ra(1)], ...
     'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
ylabel('k1');
xticklabels({'a','b','c','d','e','f','t','h','i'});
subplot(2,2,2);
bar([cost_stats_exp(2),cost_stats_mvn(2),cost_stats_mvp(2), ...
     cost_stats_m3n(2),cost_stats_m3p(2),cost_stats_m4n(2), ...
     cost_stats_m4p(2),cost_stats_rs(2),cost_stats_ra(2)], ...
     'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
ylabel('k2');
xticklabels({'a','b','c','d','e','f','t','h','i'});
subplot(2,2,3);
bar([cost_stats_exp(3),cost_stats_mvn(3),cost_stats_mvp(3), ...
     cost_stats_m3n(3),cost_stats_m3p(3),cost_stats_m4n(3), ...
     cost_stats_m4p(3),cost_stats_rs(3),cost_stats_ra(3)], ...
     'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
ylabel('k3');
xticklabels({'a','b','c','d','e','f','t','h','i'});
subplot(2,2,4);
bar([cost_stats_exp(4) cost_stats_mvn(4),cost_stats_mvp(4), ...
     cost_stats_m3n(4),cost_stats_m3p(4),cost_stats_m4n(4), ...
     cost_stats_m4p(4),cost_stats_rs(4),cost_stats_ra(4)], ...
     'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
ylabel('k4');
xticklabels({'a','b','c','d','e','f','t','h','i'});

file = '../../dhri/drafts/onILEQR/figures/simulation_stats.eps';
set(figure_cost_stas,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(figure_cost_stas, '-depsc2','-painters', file);