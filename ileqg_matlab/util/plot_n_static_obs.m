function [p_handle] = plot_n_static_obs(sobs, plot_color)

hold on;
p_handle = [];
for i = 1:length(sobs)
    for j = 1:length(sobs{i})
        p_handle = [p_handle ...
            plot(sobs{i}{j}(1), sobs{i}{j}(2), '.','color', plot_color)];
    end
end
