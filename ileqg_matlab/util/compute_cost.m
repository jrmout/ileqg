function [cost] = compute_cost(fnCost, x_traj, u_traj, dt)
%STOCHASTIC_SIMULATION 
    T = size(x_traj,2);
    cost = 0;
    for k=1:T-1
         cost = cost + dt*fnCost(x_traj(:,k), u_traj(:,k), k);
    end
    cost = cost + fnCost(x_traj(:,T), u_traj(:,T), T);
end