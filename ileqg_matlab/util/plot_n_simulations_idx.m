function [p_handle] = plot_n_simulations_idx(x_trajs, idx, add_idx, plot_color)

hold on;
p_handle = [];
for i = 1:length(x_trajs)
    if ~isempty(add_idx)
        p_handle = [p_handle ...
            plot(x_trajs{i}(idx,:) + x_trajs{i}(add_idx,:), ...
                 x_trajs{i}(idx+1,:) + x_trajs{i}(add_idx+1,:), ...
                 'color', plot_color)];
    else
        p_handle = [p_handle ...
            plot(x_trajs{i}(idx,:), x_trajs{i}(idx+1,:), 'color', plot_color)];
    end
end
