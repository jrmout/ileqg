function [avg_collisions, var_collisions] = compute_collisions(x_trajs, idx_obs, rob_size)
%COST STATS
    collisions = zeros(length(x_trajs),1);
    for i=1:length(x_trajs)
        dist = sqrt(sum(x_trajs{i}(idx_obs:idx_obs+1,:).^2,1));        
        collisions(i) = sum(dist < rob_size);
    end
    avg_collisions = mean(collisions);
    var_collisions = var(collisions);
end