function [avg_cost, var_cost] = compute_cost_stats(fnCost, x_trajs, u_trajs, dt)
%COST STATS
    costs = zeros(length(x_trajs),1);
    for i=1:length(x_trajs)
        costs(i) = compute_cost(fnCost, x_trajs{i}, u_trajs{i}, dt);
    end
    avg_cost = mean(costs);
    var_cost = var(costs);
end