function [ x_trajs, u_trajs, sobs_trajs] = stochastic_n_simulations(n_sims, dt, T, x0, x, ...
                                                        u, L, fnDyn, varargin)

sobs_means = [];
sobs_vars = [];
if nargin > 8
    sobs_means = varargin{1};
    sobs_vars = varargin{2};
end
    
x_trajs = cell(n_sims,1);
u_trajs = cell(n_sims,1);
sobs_trajs = cell(n_sims,1);
for i = 1:n_sims
    % Simulate dynamics
    [x_trajs{i}, u_trajs{i}]= stochastic_simulation(x0, x, u, L, fnDyn, dt, T);

    % Simulate static uncertain elements
    if ~isempty(sobs_vars)
        for so=1:size(sobs_vars, 3)
            sobs_trajs{i}{so} = mvnrnd(sobs_means(:,so)',sobs_vars(:,:,so),1)';
        end
    end
end
