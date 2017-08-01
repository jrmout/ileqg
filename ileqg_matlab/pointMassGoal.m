close all;

% Optimization parameters
T = 50;
dt = 0.01;
theta = 1;
x0 = [0 0 0 0 0.1 0.1 1 -0.5]';
u0 = zeros(2,1);
max_iter = 1000;
uMin = -Inf;
uMax = Inf;

% Define LQ solvers
ExpectedCost = @(lqProb, varargin)kcc_multi(1, lqProb, varargin{:});
kccOptMeanNVar = @(lqProb, varargin)kcc_multi([1 -4], lqProb, varargin{:});
kccOptMeanPVar = @(lqProb, varargin)kcc_multi([1 49], lqProb, varargin{:});
kccOptMeanNSkew = @(lqProb, varargin)kcc_multi([1 0 -8], lqProb, varargin{:});
kccOptMeanPSkew = @(lqProb, varargin)kcc_multi([1 0 343], lqProb, varargin{:});
kccOptMeanNKurt = @(lqProb, varargin)kcc_multi([1 0 0 -16], lqProb, varargin{:});
kccOptMeanPKurt = @(lqProb, varargin)kcc_multi([1 0 0 2401], lqProb, varargin{:});
kccOptMeanPVarPKur = @(lqProb, varargin)kcc_multi([1 0.0001 0.0001], lqProb, varargin{:});
kccOptMeanNVarPKur = @(lqProb, varargin)kcc_multi([1 -0.0001 0.00001], lqProb, varargin{:});
riskAverse = @(lqProb, varargin)leqr_multi([10], lqProb, varargin{:});
riskSeeking = @(lqProb, varargin)leqr_multi([-20], lqProb, varargin{:});

% Solve iteratively
[x_exp, u_exp, L_exp, cost_exp, eps_exp] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, ExpectedCost, dt, T, x0, u0, uMin, uMax, max_iter);
[x_mvn, u_mvn, L_mvn, cost_mvn, eps_mvn] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanNVar, dt, T, x0, u_exp, uMin, uMax, max_iter);
[x_mvp, u_mvp, L_mvp, cost_mvp, eps_mvp] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanPVar, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_msn, u_msn, L_msn, cost_msn, eps_msn] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanNSkew, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_msp, u_msp, L_msp, cost_msp, eps_msp] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanPSkew, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_mkn, u_mkn, L_mkn, cost_mkn, eps_mkn] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanNKurt, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_mkp, u_mkp, L_mkp, cost_mkp, eps_mkp] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, kccOptMeanPKurt, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_ra, u_ra, L_ra, cost_ra, eps_ra] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, riskAverse, dt, T, x0, u_exp, uMin, uMax, max_iter);
% [x_rs, u_rs, L_rs, cost_rs, eps_rs] = iterativeLQSolver(@pointMass_dyngoal_dyn, @pointMass_dyngoal_cost, riskSeeking, dt, T, x0, u_exp, uMin, uMax, max_iter);


%% PLOT

xmin = -0.1;
xmax = 0.65;
ymin = -0.3;
ymax = 0.15;

fig_dyngoal  = figure('units','normalized','outerposition',[0 0 0.5 0.3]);


ax1 = axes('Parent',fig_dyngoal,'Position',[0.057 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoal_plot(x_exp, u_exp, L_exp, ax1);
xlim(ax1,[xmin xmax]);
ylim(ax1,[ymin ymax]);
box(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'ax1');
legend1 = legend(ax1,'show', {'l11111111111','l2222222222', 'l3333333333'});
legend('boxoff');
set(legend1,'Position',[0.0684300341296929 0.139650872817955 0.13839590443686 0.316084788029925],'FontSize',13);


ax2 = axes('Parent',fig_dyngoal,'Position',[0.245 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
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
pointMass_dyngoal_plot(x_mvp, u_mvp, L_mvp, ax3);
xlim(ax3,[xmin xmax]);
ylim(ax3,[ymin ymax]);
box(ax3, 'on');
grid(ax3, 'on');
xlabel(ax3, 'ax3');

ax4 = axes('Parent',fig_dyngoal,'Position',[0.439 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoal_plot(x_msn, u_msn, L_msn, ax4);
set(ax4,'YTickLabel',{'','','','',''})
xlim(ax4,[xmin xmax]);
ylim(ax4,[ymin ymax]);
box(ax4, 'on');
grid(ax4, 'on');
xlabel(ax4, 'ax4');

ax5 = axes('Parent',fig_dyngoal,'Position',[0.439 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoal_plot(x_msp, u_msp, L_msp, ax5);
set(ax5,'YTickLabel',{'','','','',''})
xlim(ax5,[xmin xmax]);
ylim(ax5,[ymin ymax]);
box(ax5, 'on');
grid(ax5, 'on');
xlabel(ax5, 'ax5');

ax6 = axes('Parent',fig_dyngoal,'Position',[0.63 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoal_plot(x_mkn, u_mkn, L_mkn, ax6);
set(ax6,'YTickLabel',{'','','','',''})
xlim(ax6,[xmin xmax]);
ylim(ax6,[ymin ymax]);
box(ax6, 'on');
grid(ax6, 'on');
xlabel(ax6, 'ax6');

ax7 = axes('Parent',fig_dyngoal,'Position',[0.63 0.08 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
pointMass_dyngoal_plot(x_mkp, u_mkp, L_mkp, ax7);
set(ax7,'YTickLabel',{'','','','',''})
xlim(ax7,[xmin xmax]);
ylim(ax7,[ymin ymax]);
box(ax7, 'on');
grid(ax7, 'on');
xlabel(ax7, 'ax7');

ax8 = axes('Parent',fig_dyngoal,'Position',[0.82 0.57 0.163 0.43],...
    'PlotBoxAspectRatio',[1.66666666666667 1 4.44444444444444],...
    'DataAspectRatio',[1 1 1], 'FontSize',13);
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
pointMass_dyngoal_plot(x_ra, u_ra, L_ra, ax9);
set(ax9,'YTickLabel',{'','','','',''})
xlim(ax9,[xmin xmax]);
ylim(ax9,[ymin ymax]);
box(ax9, 'on');
grid(ax9, 'on');
xlabel(ax9, 'ax9');


file = '../drafts/onILEQR/figures/simulation_dyngoal.eps';
set(fig_dyngoal,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_dyngoal, '-depsc2','-painters', file);