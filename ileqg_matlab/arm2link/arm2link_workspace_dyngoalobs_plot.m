function [] = arm2link_workspace_dyngoalobs_plot( x, ~ )
%COST_2D_ORIENT Summary of this function goes here
%   Detailed explanation goes here

% Define Colors
black = [0 0 0];
red_averse = [219/255 70/255 70/255];
green_seeking = [60/255 169/255 60/255];
darkmagenta = [139/255 0 139/255];
white = [1 1 1];


% Arm model parameters
l1 = 0.3;
l2 = 0.33;

T = size(x,2);

plotStep = 10;
nPlots = T/plotStep;

x_goal = x(5:8,:) + x(1:4,:);
x_obs = x(9:12,:) + x(1:4,:);

clf;

for i = 1:T
    D = (x(1,i)^2+x(2,i)^2-l1^2-l2^2)/(2*l1*l2);

    if (1-D^2) > 0
        % this assumes elbow down position
        theta(2) = atan2(-sqrt(1-D^2),D);
        theta(1) = atan2(x(2,i),x(1,i)) - atan2(l2*sin(theta(2)),l1+l2*cos(theta(2)));
        
    else
        % Out of the workspace...
        theta(2) = 0;
        theta(1) = 0;
    end
    
    
    
    % Base coordinates
    taskCoord(:,1) = [0 0]';

    % First link task coordinates
    taskCoord(:,2) = [l1*cos(theta(1)) l1*sin(theta(1))];

    % EF task coordinates
    taskCoord(:,3) = [l1*cos(theta(1)) + l2*cos(theta(1) + theta(2)) l1*sin(theta(1)) + l2*sin(theta(1) + theta(2))];
    
    if mod(i,plotStep) == 0
        % Plot first link 
        pseudoGrad = ((i/plotStep)/nPlots);
        hold on;
        if (pseudoGrad == 1)
            plot(taskCoord(1,:), taskCoord(2,:), 'LineWidth', 5, 'color', black);
            plot(x_goal(1,i), x_goal(2,i), 'x', 'Markersize', 6, 'LineWidth', 5, 'color', darkmagenta);
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', red_averse);
        else
            plot(taskCoord(1,:), taskCoord(2,:), 'LineWidth', 5, 'color', pseudoGrad*green_seeking + (1 - pseudoGrad)*(white));
            plot(x_goal(1,i), x_goal(2,i), 'x', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*darkmagenta + (1 - pseudoGrad)*(white));
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*red_averse + (1 - pseudoGrad)*(white));
        end
    end
    
end

axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]) ;
drawnow;

end