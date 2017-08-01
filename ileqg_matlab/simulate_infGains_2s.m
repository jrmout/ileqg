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
Q(2,2) = 0;
R = eye(length(A));
theta_neutral = 0;
theta_averse = -0.001;
theta_seeking = 0.001;

lqProbInf =  {Q, R, A, B, theta_neutral, Gamma, Sigma};


%% 2 uncertainty sources

Gamma1 = eye(length(A));
Gamma2 = eye(length(A));
Sigma1 = eye(length(A));
Sigma2 = eye(length(A));
Gamma_all(:,:,1) = Gamma1;
Gamma_all(:,:,2) = Gamma2;
Sigma_all(:,:,1) = Sigma1;
Sigma_all(:,:,2) = Sigma2;
lqProbInf{6} = Gamma_all;
lqProbInf{7} = Sigma_all;


gamma_factor = 0.0001*(1:5:10);

L_risks = zeros(Dim*2, Dim*2, length(gamma_factor), length(gamma_factor));
thetas = [0.002 -0.002]';
lqProbInf{5} = thetas;

X_risk = [];
Y_risk = [];
Z_risk = [];
for i=1:length(gamma_factor)
    Xi = [];
    Yi = [];
    Zi = [];
    i
    for j=1:length(gamma_factor)
        Gamma_all(:,:,1) = gamma_factor(i)*Gamma1;
        Gamma_all(:,:,2) = gamma_factor(j)*Gamma2;
        lqProbInf{6} = Gamma_all;
        L_risks(:,:,i,j) = leqr_fb_multi_inf(lqProbInf);
        Xi = [Xi i];
        Yi = [Yi j];
        Zi = [Zi abs(L_risks(2,1,i,j))];
    end
    X_risk = [ X_risk ; Xi ];
    Y_risk = [ Y_risk ; Yi ];
    Z_risk = [ Z_risk ; Zi ];
end


gamma = 0.03163;

factor = 2;

gammas = [1 factor*gamma^2 ; 1 -factor*gamma^2];
lqProbInf{5} = gammas;

X_mvar = [];
Y_mvar = [];
Z_mvar = [];
for i=1:length(gamma_factor)
    Xi = [];
    Yi = [];
    Zi = [];
    i
    for j=1:length(gamma_factor)
        Gamma_all(:,:,1) = gamma_factor(i)*Gamma1;
        Gamma_all(:,:,2) = gamma_factor(j)*Gamma2;
        lqProbInf{6} = Gamma_all;
        L_risks(:,:,i,j) = kcc_fb_multi_inf(lqProbInf);
        Xi = [Xi i];
        Yi = [Yi j];
        Zi = [Zi abs(L_risks(2,1,i,j))];
    end
    X_mvar = [ X_mvar ; Xi ];
    Y_mvar = [ Y_mvar ; Yi ];
    Z_mvar = [ Z_mvar ; Zi ];
end

gammas = [1 0 factor*gamma^3 ; 1 0 -factor*gamma^3];
lqProbInf{5} = gammas;

X_mskew = [];
Y_mskew = [];
Z_mskew = [];
for i=1:length(gamma_factor)
    Xi = [];
    Yi = [];
    Zi = [];
    i
    for j=1:length(gamma_factor)
        Gamma_all(:,:,1) = gamma_factor(i)*Gamma1;
        Gamma_all(:,:,2) = gamma_factor(j)*Gamma2;
        lqProbInf{6} = Gamma_all;
        L_risks(:,:,i,j) = kcc_fb_multi_inf(lqProbInf);
        Xi = [Xi i];
        Yi = [Yi j];
        Zi = [Zi abs(L_risks(2,1,i,j))];
    end
    X_mskew = [ X_mskew ; Xi ];
    Y_mskew = [ Y_mskew ; Yi ];
    Z_mskew = [ Z_mskew ; Zi ];
end


gammas = [1 0 0 factor*gamma^4 ; 1 0 0 -factor*gamma^4];
lqProbInf{5} = gammas;

X_mkur = [];
Y_mkur = [];
Z_mkur = [];
for i=1:length(gamma_factor)
    Xi = [];
    Yi = [];
    Zi = [];
    i
    for j=1:length(gamma_factor)
        Gamma_all(:,:,1) = gamma_factor(i)*Gamma1;
        Gamma_all(:,:,2) = gamma_factor(j)*Gamma2;
        lqProbInf{6} = Gamma_all;
        L_risks(:,:,i,j) = kcc_fb_multi_inf(lqProbInf);
        Xi = [Xi i];
        Yi = [Yi j];
        Zi = [Zi abs(L_risks(2,1,i,j))];
    end
    X_mkur = [ X_mkur ; Xi ];
    Y_mkur = [ Y_mkur ; Yi ];
    Z_mkur = [ Z_mkur ; Zi ];
end

%Compute colormap
greengrayred = zeros(512,3);

% From green-seeking to gray neutral to red_averse
for i = 1:256
    pseudoGrad = (i/256);
    greengrayred(i,:) = (1 - pseudoGrad)*green_seeking + pseudoGrad*(grey_neutral);
end
for i = 1:256
    pseudoGrad = (i/256);
    greengrayred(256 + i,:) = (1 - pseudoGrad)*grey_neutral + pseudoGrad*(red_averse);
end


%% PLOT

% Expected gain is 0.998
coloraxes = [0.988 1.008];

fig  = figure('units','normalized','outerposition',[0 0 0.51 0.29]);
axes1 = axes('Parent',fig,'Position',[0.0502059202059202 0.239401496259352 0.2 0.7],'CLim',coloraxes, 'Fontsize', 14);
axes(axes1);
imagesc(flipud(Z_risk), 'CDataMapping','scaled'); 
set(axes1,'YDir','normal')
caxis( axes1, coloraxes); c=colormap(axes1, greengrayred);shading interp;
%uimage(X_risk(:,1),X_risk(:,1),Z_risk); caxis( axes1, coloraxes  );c=colormap(axes1, greengrayred);shading interp;
xlabel(axes1,'x1');
ylabel(axes1,'y');
set(axes1,'XTickLabel',{'0.2','0.4','0.8','1'}, 'Fontsize', 14)
set(axes1,'YTickLabel',{'0.2','0.4','0.8','1'}, 'Fontsize', 14)


axes2 = axes('Parent',fig,'Position',[0.279488964648537 0.239401496259352 0.2 0.7],'CLim',coloraxes, 'Fontsize', 14);
axes(axes2);
image(flipud(Z_mvar), 'CDataMapping','scaled');
set(axes2,'YDir','normal')
caxis( axes2, coloraxes); c=colormap(axes2, greengrayred);shading interp;
%pcolor(axes2, X_mvar,Y_mvar,Z_mvar); caxis( coloraxes  );colormap(c); shading interp;
xlabel(axes2,'x2');
set(axes2,'YTickLabel',{'','','','',''})
set(axes2,'XTickLabel',{'0.2','0.4','0.8','1'}, 'Fontsize', 14)


axes3 = axes('Parent',fig,'Position',[0.509652509652508 0.239401496259352 0.2 0.7],'CLim',coloraxes, 'Fontsize', 14);
axes(axes3);
image(flipud(Z_mskew), 'CDataMapping','scaled'); 
set(axes3,'YDir','normal')
caxis( axes3, coloraxes); c=colormap(axes3, greengrayred);shading interp;
%pcolor(axes3, X_mskew,Y_mskew,Z_mskew); caxis( coloraxes  );colormap(c); shading interp;
xlabel(axes3,'x3');
set(axes3,'YTickLabel',{'','','','',''})
set(axes3,'XTickLabel',{'0.2','0.4','0.8','1'}, 'Fontsize', 14)


axes4 = axes('Parent',fig,'Position',[0.74002574002574 0.239401496259352 0.2 0.7],'CLim',coloraxes, 'Fontsize', 14);
axes(axes4);
image(flipud(Z_mkur), 'CDataMapping','scaled'); 
set(axes4,'YDir','normal')
caxis( axes4, coloraxes); c=colormap(axes4, greengrayred);shading interp;
%pcolor(axes4, X_mkur,Y_mkur,Z_mkur); caxis( coloraxes  ); colormap(c); shading interp;
xlabel(axes4,'x4');
set(axes4,'YTickLabel',{'','','','',''})
set(axes4,'XTickLabel',{'0.2','0.4','0.8','1'}, 'Fontsize', 14)
colorbar(axes4, 'Position', [0.958655083655083 0.234413965087282 0.0160875160875161 0.703241895261845], 'Fontsize', 14);
set(axes4,'Position',[0.74002574002574 0.239401496259352 0.2 0.7]);


% Create textbox
annotation(fig,'textbox',[0.0380180180180166 0.163 0.0186615186615195 0.0688254364089775],'String',{'0'},'FontSize',14,'FitBoxToText','off','LineStyle','none');
annotation(fig,'textbox',[0.272947232947227 0.163 0.0186615186615195 0.0688254364089775],'String',{'0'},'FontSize',14,'FitBoxToText','off','LineStyle','none');
annotation(fig,'textbox',[0.734440154440125 0.163 0.0186615186615195 0.0688254364089775],'String',{'0'},'FontSize',14,'FitBoxToText','off','LineStyle','none');
annotation(fig,'textbox',[0.50337194337193 0.163 0.0186615186615195 0.0688254364089775],'String',{'0'},'FontSize',14,'FitBoxToText','off','LineStyle','none');
annotation(fig,'textbox',[0.140283140283139 0.02 0.0186615186615195 0.0688254364089775],'String',{'p1'},'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
annotation(fig,'textbox',[0.369420849420848 0.02 0.0186615186615195 0.0688254364089775],'String',{'p2'},'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
annotation(fig,'textbox',[0.603655083655082 0.02 0.0186615186615195 0.0688254364089775],'String',{'p3'},'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
annotation(fig,'textbox',[0.834028314028313 0.02 0.0186615186615195 0.0688254364089775],'String',{'p4'},'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
annotation(fig,'textbox',[0.960797940797937 0.05 0.0186615186615195 0.0688254364089775],'String',{'g'},'FontSize',14,'FitBoxToText','off','LineStyle','none');


file = '../drafts/onILEQR/figures/simulation_infGains_2s.eps';
filepdf = '../drafts/onILEQR/figures/simulation_infGains_2s.pdf';
set(fig,'PaperPositionMode','auto');
iptsetpref('ImshowBorder','tight');
%export_fig '../drafts/onILEQR/figures/simulation_infGains_2s.pdf';
print(fig, '-depsc2','-painters', file);
%fix_pcolor_eps(file);

