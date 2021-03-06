function [ ] = pointMass_plot( x, ~ , x_tar)
% Plots a point mass

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

% x_tar is the goal in task space

plotStep = 10;
nPlots = T/plotStep;

clf;

for i = 1:T
    
    if mod(i,plotStep) == 0
        pseudoGrad = ((i/plotStep)/nPlots);
        hold on;
        if (pseudoGrad == 1)
            plot(x(1,i), x(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', black);
        else
            plot(x(1,i), x(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*green_seeking + (1 - pseudoGrad)*(white));
        end
    end
    
end

% Plot target
plot(x_tar(1), x_tar(2), '+', 'LineWidth', 5, 'color', darkred);
axis([-1 1 -1 1]) ;

drawnow;

end

