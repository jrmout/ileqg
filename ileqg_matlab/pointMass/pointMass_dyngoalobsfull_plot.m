function [ ] = pointMass_dyngoalobsfull_plot( x, ~, L, sobs_pos, sobs_Sigma, varargin)
% Plots a point mass

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

% x_tar is the goal in task space

plotStep = 2;
nPlots = T/plotStep;



x_rob = x(1:4,:);
x_goal = x(5:8,:) + x(1:4,:);
x_obs = x(9:12,:) + x(1:4,:);

L_rob = L(:,1:4,:);
L_goal = L(:,5:8,:);
L_obs = L(:,9:12,:);


hold on;

% Plot static Obstacles
for i = 1:size(sobs_pos,2)
    [V,D] = eig(sobs_Sigma(:,:,i));    
    %Scale eigenvectors
    Di = diag(D);
    Di = Di / 10;
    if ~(any(~isreal(Di))) && ~any(isinf(Di)) && ~any(isnan(Di))
        static_obstacle = plot(sobs_pos(1,i), sobs_pos(2,i), 'o', 'Markersize', 2, 'LineWidth', 5,'color', color_static_obs, 'Parent', axes_handle);
        set(static_obstacle, 'Tag','static obstacle');
        %[el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), V(:,1)*Di(1), V(:,2)*Di(2), 50);
        %plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color',0.8*color_static_obs + 0.2*(white), 'Parent', axes_handle);                
        %[el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), 2*V(:,1)*Di(1), 2*V(:,2)*Di(2), 50);
        %plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color', 0.6*color_static_obs + 0.4*(white), 'Parent', axes_handle);                
        %[el1_goal,el2_goal] = computeEllipse(sobs_pos(1,i), sobs_pos(2,i), 3*V(:,1)*Di(1), 3*V(:,2)*Di(2), 50);
        %plot(el1_goal, el2_goal, '-', 'LineWidth', 2, 'color', 0.3*color_static_obs + 0.7*(white), 'Parent', axes_handle);                
        
    end
end

for i = 1:T
    
    if mod(i,plotStep) == 0
        pseudoGrad = ((i/plotStep)/nPlots);
        if (i == plotStep)
            plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', black, 'Parent', axes_handle);
            plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', black, 'Parent', axes_handle);
            plot(x_rob(1,i), x_rob(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', black, 'Parent', axes_handle);
        else
            if (pseudoGrad < 1)
                [V,D] = eig(L_obs(:,1:2,i));
                Di = diag(D);
                %Scale eigenvectors
                Di = Di / 1000;
                if ~(any(~isreal(Di))) && ~any(isinf(Di)) && ~any(isnan(Di))
                    %[el1_obs,el2_obs] = computeEllipse(x_rob(1,i), x_rob(2,i), V(:,1)*Di(1), V(:,2)*Di(2), 50);
                    %gain_obstacle = plot(el1_obs, el2_obs, '-', 'LineWidth', 1, 'color', pseudoGrad*color_obstacle_gain + (1 - pseudoGrad)*(light_color_obstacle_gain), 'Parent', axes_handle);            
                    %set(gain_obstacle, 'Tag','gain obstacle');
                end

                [V,D] = eig(L_goal(:,1:2,i));
                %Scale eigenvectors
                Di = diag(D);
                Di = Di / 1000;           
                if ~(any(~isreal(Di))) && ~any(isinf(Di)) && ~any(isnan(Di))
                    %[el1_goal,el2_goal] = computeEllipse(x_rob(1,i), x_rob(2,i), V(:,1)*Di(1), V(:,2)*Di(2), 50);
                    %plot(el1_goal, el2_goal, '-', 'LineWidth', 1, 'color', pseudoGrad*color_goal_gain + (1 - pseudoGrad)*(light_color_goal_gain), 'Tag','gain goal', 'Parent', axes_handle);                
                end
            end
            
            goal = plot(x_goal(1,i), x_goal(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', pseudoGrad*color_goal + (1 - pseudoGrad)*(light_color_goal), 'Parent', axes_handle);
            set(goal, 'Tag','goal');
            dynobs = plot(x_obs(1,i), x_obs(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', pseudoGrad*color_obstacle + (1 - pseudoGrad)*(light_color_obstacle), 'Parent', axes_handle);
            set(dynobs, 'Tag','dynamic obstacle');
            robot = plot(x_rob(1,i), x_rob(2,i), 'o', 'Markersize', 2, 'LineWidth', 5, 'color', pseudoGrad*color_robot + (1 - pseudoGrad)*(light_color_robot), 'Parent', axes_handle);
            set(robot, 'Tag','robot');


        end
    end
end

drawnow;

end

