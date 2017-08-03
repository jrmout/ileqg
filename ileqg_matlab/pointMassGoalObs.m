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
x_goal = [0.1 0.1 1 -0.5]';
x_obs = [0.15 0.15 1 -0.5]';
x0 = [x_rob ; (x_goal - x_rob) ; (x_obs - x_rob)];

gamma_neg = 4;
gamma_pos = 5;

% Dynamics
dynamics_function = @pointMass_dyngoalobs_dyn;

% Cost
cost_function = @pointMass_dyngoalobs_cost;

% Plot (during optimization)
plot_function = @pointMass_dyngoalobs_plot;


% LQ solvers
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
kccOptMeanNVar = @(lqProb, varargin)kcc_multi(...
                                [1 -(gamma_neg^2)/2], lqProb, varargin{:});
kccOptMeanPVar = @(lqProb, varargin)kcc_multi(...
                                [1 gamma_pos^2 /2], lqProb, varargin{:});
kccOptMeanNThird = @(lqProb, varargin)kcc_multi(...
                                [1 0 -(gamma_neg^3) / 6], lqProb, varargin{:});
kccOptMeanPThird = @(lqProb, varargin)kcc_multi(...
                                [1 0 gamma_pos^3 / 6], lqProb, varargin{:});
kccOptMeanNFourth = @(lqProb, varargin)kcc_multi(...
                                [1 0 0 -(gamma_neg^4) / 24], lqProb, varargin{:});
kccOptMeanPFourth = @(lqProb, varargin)kcc_multi(...
                                [1 0 0 gamma_pos^4 / 24], lqProb, varargin{:});
riskAverse = @(lqProb, varargin)leqr_multi([4], lqProb, varargin{:});
riskSeeking = @(lqProb, varargin)leqr_multi([-50], lqProb, varargin{:});

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
 [x_msn, u_msn, L_msn, cost_msn, eps_msn] = ...
     iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNThird, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
 [x_msp, u_msp, L_msp, cost_msp, eps_msp] = ...
     iterativeLQSolver(dynamics_function, cost_function, kccOptMeanPThird, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
 [x_mkn, u_mkn, L_mkn, cost_mkn, eps_mkn] = ...
     iterativeLQSolver(dynamics_function, cost_function, kccOptMeanNFourth, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
 [x_mkp, u_mkp, L_mkp, cost_mkp, eps_mkp] = ...
     iterativeLQSolver(dynamics_function, cost_function, kccOptMeanPFourth, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
 [x_ra, u_ra, L_ra, cost_ra, eps_ra] = ...
     iterativeLQSolver(dynamics_function, cost_function, riskAverse, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);
 [x_rs, u_rs, L_rs, cost_rs, eps_rs] = ...
     iterativeLQSolver(dynamics_function, cost_function, riskSeeking, ...
                       dt, T, x0, u_exp, uMin, uMax, max_iter, plot_function);

%% Simulate solutions and get statistics

                   
%% PLOT
xmin = -0.15;
xmax = 0.6;
ymin = -0.275;
ymax = 0.175;

fig_dyngoal  = figure('units','normalized','outerposition',[0 0 0.5 0.3]);


ax1 = axes('Parent',fig_dyngoal,'Position',[0.057 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_exp, u_exp, L_exp, ax1);
xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
legend1 = legend(ax1,'show', {'l11111111111','l2222222222', 'l3333333333', 'l444444444', 'l55555555'});
legend('boxoff');
set(legend1,'Position',[0.0684300341296929 0.139650872817955 0.13839590443686 0.316084788029925],'FontSize',13);


ax2 = axes('Parent',fig_dyngoal,'Position',[0.245 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_mvn, u_mvn, L_mvn, ax2);
xlim(ax2,[xmin xmax]);
ylim(ax2,[ymin ymax]);
set(ax2,'YTickLabel',{'','','','',''})
box(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'ax2');

ax3 = axes('Parent',fig_dyngoal,'Position',[0.245 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_mvp, u_mvp, L_mvp, ax3);
xlim(ax3,[xmin xmax]);
ylim(ax3,[ymin ymax]);
box(ax3, 'on');
grid(ax3, 'on');
xlabel(ax3, 'ax3');

ax4 = axes('Parent',fig_dyngoal,'Position',[0.439 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_msn, u_msn, L_msn, ax4);
set(ax4,'YTickLabel',{'','','','',''})
xlim(ax4,[xmin xmax]);
ylim(ax4,[ymin ymax]);
box(ax4, 'on');
grid(ax4, 'on');
xlabel(ax4, 'ax4');

ax5 = axes('Parent',fig_dyngoal,'Position',[0.439 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_msp, u_msp, L_msp, ax5);
set(ax5,'YTickLabel',{'','','','',''})
xlim(ax5,[xmin xmax]);
ylim(ax5,[ymin ymax]);
box(ax5, 'on');
grid(ax5, 'on');
xlabel(ax5, 'ax5');

ax6 = axes('Parent',fig_dyngoal,'Position',[0.63 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_mkn, u_mkn, L_mkn, ax6);
set(ax6,'YTickLabel',{'','','','',''})
xlim(ax6,[xmin xmax]);
ylim(ax6,[ymin ymax]);
box(ax6, 'on');
grid(ax6, 'on');
xlabel(ax6, 'ax6');

ax7 = axes('Parent',fig_dyngoal,'Position',[0.63 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_mkp, u_mkp, L_mkp, ax7);
set(ax7,'YTickLabel',{'','','','',''})
xlim(ax7,[xmin xmax]);
ylim(ax7,[ymin ymax]);
box(ax7, 'on');
grid(ax7, 'on');
xlabel(ax7, 'ax7');

ax8 = axes('Parent',fig_dyngoal,'Position',[0.82 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_rs, u_rs, L_rs, ax8);
set(ax8,'YTickLabel',{'','','','',''})
xlim(ax8,[xmin xmax]);
ylim(ax8,[ymin ymax]);
box(ax8, 'on');
grid(ax8, 'on');
xlabel(ax8, 'ax8');

ax9 = axes('Parent',fig_dyngoal,'Position',[0.82 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoalobs_plot(x_ra, u_ra, L_ra, ax9);
set(ax9,'YTickLabel',{'','','','',''})
xlim(ax9,[xmin xmax]);
ylim(ax9,[ymin ymax]);
box(ax9, 'on');
grid(ax9, 'on');
xlabel(ax9, 'ax9');


file = '../drafts/onILEQR/figures/simulation_dyngoalobs.eps';
set(fig_dyngoal,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_dyngoal, '-depsc2','-painters', file);