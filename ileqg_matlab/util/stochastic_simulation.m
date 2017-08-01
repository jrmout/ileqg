function [ x_traj, u_traj] = stochastic_simulation( x0, x_des, l, L, fnDyn, dt, T)
%STOCHASTIC_SIMULATION 
    x = x0;
    x_traj = zeros(size(x,1),T);
    u_traj = zeros(size(l,1),T);
    x_traj(:,1) = x;
    for k=1:T-1
        u = l(:,k) - L(:,:,k)*(x_des(:,k) - x);
        [x_dot, ~, ~, Gamma_k, Sigma_k] = fnDyn(x, u);
        noise = zeros(size(x0));
        for s=1:size(Gamma_k,3)
            noise = noise + mvnrnd(zeros(size(x)),Gamma_k(:,:,s)*Sigma_k(:,:,s),1)';
        end
        x = x + dt*(x_dot + noise);
        x_traj(:,k+1) = x;
        u_traj(:,k) = u;
    end
end