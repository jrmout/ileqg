function [sample_cumulants] = compute_cost_stats(fnCost, x_trajs, u_trajs, dt)
%COST STATS
    costs = zeros(length(x_trajs),1);
    for i=1:length(x_trajs)
        costs(i) = compute_cost(fnCost, x_trajs{i}, u_trajs{i}, dt);
    end
    % Sample first cumulant (mean)
    sample_cumulants(1) = mean(costs);
    
    % Sample second cumulant (variance)
    sample_cumulants(2) = sum((costs - sample_cumulants(1)).^2)/size(costs,1);
    
    % Sample third cumulant
    sample_cumulants(3) = sum((costs - sample_cumulants(1)).^3)/size(costs,1);
    
    % Sample fourth cumulant
    sample_cumulants(4) = sum((costs - sample_cumulants(1)).^4)/size(costs,1) ...
                            - 3 * sample_cumulants(2)^2;
end