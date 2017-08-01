function [ ] = carLikeRobot_dyngoalobs_plot( x, ~ )
% Plots a point carLike robot

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


T = size(x,2);

% x_tar is the goal state

plotStep = 10;
nPlots = T/plotStep;

clf;

arrowLength = 0.1;

x_goal = x(6:10,:) + x(1:5,:);
x_obs = x(11:12,:) + x(1:2,:);

for i = 1:T
    
    if mod(i,plotStep) == 0
        % Plot first link 
        pseudoGrad = ((i/plotStep)/nPlots);
        hold on;
        if (pseudoGrad == 1)
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', red_averse);
            
            quiver(x_goal(1,i), x_goal(2,i), arrowLength*cos(x_goal(3,i)), arrowLength*sin(x_goal(3,i)), '-','LineWidth',4, 'Color', darkmagenta, 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3);
            plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', darkmagenta);
            
            quiver(x(1,i), x(2,i), arrowLength*cos(x(3,i)), arrowLength*sin(x(3,i)), '-','LineWidth',4, 'Color', black, 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3);
            plot(x(1,i), x(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', black);
        else
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*red_averse + (1 - pseudoGrad)*(white));
            
            quiver(x_goal(1,i), x_goal(2,i), arrowLength*cos(x_goal(3,i)), arrowLength*sin(x_goal(3,i)), '-','LineWidth',4, 'Color', pseudoGrad*darkgrey + (1 - pseudoGrad)*(white), 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3);
            plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*darkmagenta + (1 - pseudoGrad)*(white));
            
            quiver(x(1,i), x(2,i), arrowLength*cos(x(3,i)), arrowLength*sin(x(3,i)), '-','LineWidth',4, 'Color', pseudoGrad*orange + (1 - pseudoGrad)*(white), 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3);
            plot(x(1,i), x(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*green_seeking + (1 - pseudoGrad)*(white));
        end
    end
    
end
axis([-1 1 -1 1]) ;
drawnow;

end

