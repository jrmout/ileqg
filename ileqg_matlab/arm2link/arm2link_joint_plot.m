function [] = arm2link_joint_plot( x, u, x_tar)
%COST_2D_ORIENT Summary of this function goes here
%   Detailed explanation goes here

% My favorite colors 
black = [0 0 0];
red_averse = [219/255 70/255 70/255];
grey_neutral = [0.4 0.4 0.4];
darkgrey = [0.3 0.3 0.3];
green_seeking = [60/255 169/255 60/255];
darkblue = [0/255 120/255 250/255];
orange = [0.9 120/255 0];
darkmagenta = [139/255 0 139/255];
darkgreen = green_seeking;
darkred = red_averse;
white = [1 1 1];

% Model parameters
m1 = 1.4;
m2 = 1.1;
l1 = 0.3;
l2 = 0.33;
s1 = 0.11;
s2 = 0.16;
i1 = 0.025;
i2 = 0.045;
d1 = i1 + i2 + m2*l1^2;
d2 = m2 * l1 * s2;
d3 = i2;
b11 = 0.05;
b22 = 0.05;
b12 = 0.025;
b21 = 0.025;

T = size(x,2);

% Plot arm trajectory
% Base coordinates
targetCoord(:,1) = [0 0]';

% First link task coordinates
targetCoord(:,2) = [l1*cos(x_tar(1)) l1*sin(x_tar(1))];

% EF task coordinates
targetCoord(:,3) = [l1*cos(x_tar(1)) + l2*cos(x_tar(1) + x_tar(2)) l1*sin(x_tar(1)) + l2*sin(x_tar(1) + x_tar(2))];


plotStep = 10;
nPlots = T/plotStep;

clf;

for i = 1:T
    
    % Base coordinates
    taskCoord(:,1) = [0 0]';

    % First link task coordinates
    taskCoord(:,2) = [l1*cos(x(1,i)) l1*sin(x(1,i))];

    % EF task coordinates
    taskCoord(:,3) = [l1*cos(x(1,i)) + l2*cos(x(1,i) + x(2,i)) l1*sin(x(1,i)) + l2*sin(x(1,i) + x(2,i))];
    
    if mod(i,plotStep) == 0
        % Plot first link 
        pseudoGrad = ((i/plotStep)/nPlots);
        hold on;
        if (pseudoGrad == 1)
            plot(taskCoord(1,:), taskCoord(2,:), 'LineWidth', 5, 'color', black);
        else
            plot(taskCoord(1,:), taskCoord(2,:), 'LineWidth', 5, 'color', pseudoGrad*green_seeking + (1 - pseudoGrad)*(white));
        end
    end
    
end

% Plot target
plot(targetCoord(1,:), targetCoord(2,:), 'LineWidth', 5, 'color', darkred);
axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]) ;

drawnow;

end