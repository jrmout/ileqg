function [ ] = carLikeRobot_dyngoalobsfull_plot( x, ~ , L, sobs_pos, sobs_Sigma, varargin)
% Plots a point carLike robot

% My favorite colors 
black = [0 0 0];
red_averse = [219/255 70/255 70/255];
grey_neutral = [0.4 0.4 0.4];
darkgrey = [0.3 0.3 0.3];
green_seeking = [60/255 169/255 60/255];
darkblue = [0/255 120/255 250/255];
darkmagenta = [139/255 0 139/255];
darkgreen = green_seeking;
darkred = red_averse;
white = [1 1 1];
color_robot = green_seeking;
carLength = 0.1;
carWidth = 0.05;

% Colors
color_obstacle =  [251,154,153]./255;
color_robot = [178,223,138]./255;
color_goal =  [166,206,227]./255;
color_goal_gain = [31,120,180]./255;
color_obstacle_gain = [227,26,28]./255;
color_static_obs =  [253,191,111]./255;
color_obstacle = color_obstacle*0.6 + color_obstacle_gain*0.4;
color_robot = color_robot*0.6 + [51,160,44]./255*0.4;
color_goal = color_goal* 0.6 + color_goal_gain* 0.4;
color_static_obs = color_static_obs* 0.6 + [227,26,28]./255 * 0.4;

white = [1 1 1];
black = [0.3 0.3 0.3];
light_color_goal_gain = color_goal_gain*0.9 + white*0.1;
light_color_goal = color_goal*0.9 + white*0.1;
light_color_obstacle = color_obstacle*0.9 + white*0.1;
light_color_obstacle_gain = color_obstacle_gain*0.9 + white*0.1;
light_color_robot = color_robot*0.9 + white*0.1;


if (isempty(varargin))
    clf;
    axes_handle = gca;
else
    axes_handle = varargin{1};
end

T = size(x,2);

plotStep = 40;
nPlots = T/plotStep;

arrowLength = 0.1;

x_rob = x(1:5,:);
x_goal = x(6:10,:) + x(1:5,:);
x_obs = x(11:12,:) + x(1:2,:);


L_rob = L(:,1:4,:);
L_goal = L(:,5:10,:);
L_obs = L(:,11:14,:);


hold on;

% Plot static Obstacles
for i = 1:size(sobs_pos,2)
    [V,D] = eig(sobs_Sigma(:,:,i));    
    %Scale eigenvectors
    Di = diag(D);
    Di = Di / 10;
    if ~(any(~isreal(Di))) && ~any(isinf(Di)) && ~any(isnan(Di))
        static_obstacle = plot(sobs_pos(1,i), sobs_pos(2,i), 'o', 'Markersize', 6, 'LineWidth', 5,'color', color_static_obs, 'Parent', axes_handle);
        set(static_obstacle, 'Tag','static obstacle');
        [el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), V(:,1)*Di(1), V(:,2)*Di(2), 50);
        plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color',0.8*color_static_obs + 0.2*(white), 'Parent', axes_handle);                
        [el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), 2*V(:,1)*Di(1), 2*V(:,2)*Di(2), 50);
        plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color', 0.6*color_static_obs + 0.4*(white), 'Parent', axes_handle);                
        [el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), 3*V(:,1)*Di(1), 3*V(:,2)*Di(2), 50);
        plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color', 0.3*color_static_obs + 0.7*(white), 'Parent', axes_handle);                
        
    end
end


for i = 1:T
    if (i == 1)
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', black, 'Parent', axes_handle);
            
            %quiver(x_goal(1,i), x_goal(2,i), arrowLength*cos(x_goal(3,i)), arrowLength*sin(x_goal(3,i)), '-','LineWidth',1, 'Color', black, 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3, 'Parent', axes_handle);
            %plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', black, 'Parent', axes_handle);
            
            quiver(x(1,i), x(2,i), arrowLength*cos(x(3,i)), arrowLength*sin(x(3,i)), '-','LineWidth',4, 'Color', black, 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3, 'Parent', axes_handle);
            DrawRectangle([x(1,i) x(2,i) carLength/2 carWidth-0.01 x(3,i)],'color',black, 'LineWidth', 3); 
            DrawRectangle([x(1,i) x(2,i) carLength carWidth x(3,i)],'color',black, 'LineWidth', 3);
    end
    
    if mod(i,plotStep) == 0
        % Plot first link 
        pseudoGrad = ((i/plotStep)/nPlots);
        hold on;
        plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*color_obstacle + (1 - pseudoGrad)*(light_color_obstacle), 'Parent', axes_handle, 'Tag', 'dynamic obstacle', 'Parent', axes_handle);

        quiver(x_goal(1,i), x_goal(2,i), arrowLength*cos(x_goal(3,i)), arrowLength*sin(x_goal(3,i)), '-','LineWidth',1, 'Color', pseudoGrad*color_goal + (1 - pseudoGrad)*(light_color_goal), 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3, 'Parent', axes_handle);
        plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 6, 'LineWidth', 5, 'color', pseudoGrad*color_goal + (1 - pseudoGrad)*(light_color_goal), 'Tag', 'goal', 'Parent', axes_handle);

        % Draw car
        quiver(x(1,i), x(2,i), arrowLength*cos(x(3,i)), arrowLength*sin(x(3,i)), '-','LineWidth',1, 'Color', pseudoGrad*color_robot + (1 - pseudoGrad)*(light_color_robot), 'AutoScale', 'off', 'MaxHeadSize', 0.9,'LineWidth',3, 'Parent', axes_handle);

        DrawRectangle([x(1,i) x(2,i) carLength/2 carWidth-0.01 x(3,i)],'color',pseudoGrad*color_robot + (1 - pseudoGrad)*(light_color_robot), 'LineWidth', 3); 
        DrawRectangle([x(1,i) x(2,i) carLength carWidth x(3,i)],'color',pseudoGrad*color_robot + (1 - pseudoGrad)*(light_color_robot), 'LineWidth', 3, 'Tag', 'robot');
    end
    
end
axis([-1 1 -1 1]) ;
drawnow;

end

