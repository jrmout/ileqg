%% Simulate a point mass
addpath('./lqsolutions');


% My favorite colors 
black = [0 0 0];
white = [1 1 1];
red_averse = [219/255 70/255 70/255];
grey_neutral = [0.4 0.4 0.4];
darkgrey = [0.3 0.3 0.3];
green_seeking = [60/255 169/255 60/255];
darkblue = [0/255 120/255 250/255];
orange = [0.9 120/255 0];
darkmagenta = [139/255 0 139/255];
darkgreen = green_seeking;
darkred = red_averse;


TA = 1e-3 * 5;
Dim = 1;

%% Mass and Damping matrix for the admittance
M = 1;
D = 1;

%% Discrete System Dynamics with sample time of TA (x_{k+1} = A*x_k + B*u)
A = [eye(Dim) eye(Dim)*TA;zeros(Dim) eye(Dim)-(M\D*TA)];
B = [zeros(Dim) zeros(Dim); zeros(Dim) inv(M)*TA];
Gamma = eye(length(A));
Sigma = eye(length(A));

Q = 1*eye(length(A));
Q(2,2) = 1;
R = eye(length(A));
theta_neutral = 0;
theta_averse = -0.001;
theta_seeking = 0.001;

lqProbInf =  {Q, R, A, B, theta_neutral, Gamma, Sigma};
gamma_factor = 0.00013*(1:1:100);


%% 1 single uncertainty source
L_neutral = zeros(Dim*2, Dim*2, length(gamma_factor));
L_n1 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_n2 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_n3 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_p1 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_p2 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_p3 = zeros(Dim*2, Dim*2, length(gamma_factor));


% Risk-sensitive solutions

lqProbInf{5} = 0;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_neutral(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end
disp('negative risks')
lqProbInf{5} = -0.005;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_n1(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = -0.01;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_n2(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = -0.02;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_n3(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end
disp('positive risks')
lqProbInf{5} = 0.005;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_p1(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = 0.01;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_p2(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = 0.02;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_p3(:,:,i) = leqr_fb_multi_inf(lqProbInf);
end

% Cost-cumulant solutions
L_cn1 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_cn2 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_cn3 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_cp1 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_cp2 = zeros(Dim*2, Dim*2, length(gamma_factor));
L_cp3 = zeros(Dim*2, Dim*2, length(gamma_factor));


gamma = 0.025;

lqProbInf{5} = 1;
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_neutral(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end
disp('negative cumulants')
lqProbInf{5} = [1 -gamma^2];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cn1(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = [1 0 -(gamma^3)];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cn2(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = [1 0 0 -(gamma^4)];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cn3(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end
disp('positive cumulants')
lqProbInf{5} = [1 (gamma^2)];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cp1(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = [1 0 (gamma^3)];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cp2(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end

lqProbInf{5} = [1 0 0 (gamma^4)];
for i=1:length(gamma_factor)
    lqProbInf{6} = gamma_factor(i).*Gamma;
    L_cp3(:,:,i) = kcc_fb_multi_inf(lqProbInf);
end



% Colors
green_seeking_n1 = green_seeking*0.2 + grey_neutral*0.8;
green_seeking_n2 = green_seeking*0.6 + grey_neutral*0.4;
green_seeking_n3 = green_seeking*1;

red_averse_n1 = red_averse*0.2 + grey_neutral*0.8;
red_averse_n2 = red_averse*0.6 + grey_neutral*0.4;
red_averse_n3 = red_averse*1;

fig_risks  = figure('units','normalized','outerposition',[0 0 0.2 0.3]);
axes1 = axes('Parent',fig_risks,'Position',[0.12 0.17 0.85 0.8], 'Fontsize', 16);
hold(axes1, 'on');
grid(axes1, 'on');
box(axes1, 'on');
ylim(axes1,[0.9 1.1]);
xlabel(axes1, 'x');
ylabel(axes1, 'y');
set(axes1,'XTickLabel',{'0','2','4','8','10', '12'});
plot(abs(reshape(L_neutral(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', grey_neutral);
plot(abs(reshape(L_n1(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', green_seeking_n1);
plot(abs(reshape(L_n2(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', green_seeking_n2);
plot(abs(reshape(L_n3(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', green_seeking_n3);
plot(abs(reshape(L_p1(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', red_averse_n1);
plot(abs(reshape(L_p2(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', red_averse_n2);
plot(abs(reshape(L_p3(2,1,:),[length(gamma_factor) 1])), 'Parent', axes1, 'Linewidth', 3, 'color', red_averse_n3);
annotation(fig_risks,'textbox',[0.511289198606271 0.905162218831032 0.0400696864111498 0.0538628428927678],'String',{'r3'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.712473867595818 0.905124750313272 0.0400696864111498 0.0538628428927678],'String',{'r2'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.91637630662021 0.857619632510762 0.0400696864111498 0.0538628428927678],'String',{'r1'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.518118466898954 0.17932655922383 0.0400696864111498 0.0538628428927678],'String',{'g3'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.771637630662021 0.177186635401544 0.0400696864111498 0.0538628428927678],'String',{'g2'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.924111498257839 0.259735487152763 0.0400696864111498 0.0538628428927678],'String',{'g1'},'FitBoxToText','off','LineStyle','none');
annotation(fig_risks,'textbox',[0.931846689895469 0.58079775685165 0.0400696864111498 0.0538628428927678],'String',{'e'},'FitBoxToText','off','LineStyle','none');


% Colors
green_seeking_n1 = green_seeking*0.5 + white*0.5;
green_seeking_n2 = green_seeking*0.75 + white*0.25;
green_seeking_n3 = green_seeking*1;

red_averse_n1 = red_averse*0.5 + white*0.5;
red_averse_n2 = red_averse*0.75 + white*0.25;
red_averse_n3 = red_averse*1;

fig_cumulants  = figure('units','normalized','outerposition',[0 0 0.2 0.3]);
axes2 = axes('Parent',fig_cumulants,'Position',[0.12 0.17 0.85 0.8], 'Fontsize', 16);
hold(axes2, 'on');
grid(axes2, 'on');
box(axes2, 'on');
ylim(axes2,[0.9 1.1]);
xlabel(axes2, 'x');
ylabel(axes2, 'y');
set(axes2,'XTickLabel',{'0','2','4','8','10', '12'});
plot(abs(reshape(L_neutral(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2, 'Linewidth', 3, 'color', grey_neutral);
plot(abs(reshape(L_cn1(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', green_seeking_n1);
plot(abs(reshape(L_cn2(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', green_seeking_n2);
plot(abs(reshape(L_cn3(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', green_seeking_n3);
plot(abs(reshape(L_cp1(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', red_averse_n1);
plot(abs(reshape(L_cp2(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', red_averse_n2);
plot(abs(reshape(L_cp3(2,1,:),[length(gamma_factor) 1])), 'Parent', axes2,'Linewidth', 3, 'color', red_averse_n3);
annotation(fig_cumulants,'textarrow',[0.768292682926829 0.926829268292682],[0.776119402985075 0.654228855721393],'TextEdgeColor','none','String',{'msp'});
annotation(fig_cumulants,'textarrow',[0.853658536585366 0.940766550522648],[0.917910447761194 0.895522388059702],'TextEdgeColor','none','String',{'mkp'});
annotation(fig_cumulants,'textarrow',[0.757839721254355 0.897212543554007],[0.218905472636816 0.213930348258706],'TextEdgeColor','none','String',{'mkn'});
annotation(fig_cumulants,'textarrow',[0.681184668989547 0.958188153310105],[0.325870646766169 0.450248756218905],'TextEdgeColor','none','String',{'msn'});
annotation(fig_cumulants,'textbox',[0.935481411091168 0.531796044796979 0.0373519163763067 0.056356608478803],'String',{'e'},'FitBoxToText','off','LineStyle','none');
annotation(fig_cumulants,'textarrow',[0.543554006968641 0.878048780487805],[0.460199004975124 0.517412935323383],'TextEdgeColor','none','String',{'mvn'});
annotation(fig_cumulants,'textarrow',[0.548780487804878 0.876306620209059],[0.649253731343284 0.601990049751244],'TextEdgeColor','none','String',{'mvp'});



file = '../drafts/onILEQR/figures/simulation_infGains_Risk.eps';
set(fig_risks,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_risks, '-depsc2','-painters', file);


file = '../drafts/onILEQR/figures/simulation_infGains_Cumulants.eps';
set(fig_cumulants,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
print(fig_cumulants, '-depsc2','-painters', file);