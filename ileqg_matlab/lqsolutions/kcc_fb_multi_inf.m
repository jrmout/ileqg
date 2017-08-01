function [L] = kcc_fb_multi_inf(lqProbInf, varargin)

% Approximates optimal feedback gain numerically for an infinite horizon

% -- INPUTS
% LqProbInf is a cell with:
% - Q                     Quadratic State Weighting matrix - n x n 
% - R                     Quadratic Control input weighting - m x m 
% - A                     System matrix - n x n x T , where n is the dimension of the state
% - B                     Control matrix - n x m x T , where m is the control input dimension
% - theta                 Risk-sensitivities S
% - Gamma                 Control noise matrix matrix elements - n x qs x S noise inputs
% - Sigma                 covariance matrix - qs x qs x S

% theta                 risk sensitivities S

% -- OUTPUTS
% L                     feedback control component  - m x n x T-1


%% INITIALIZATION
Q = lqProbInf{1};
R = lqProbInf{2};
A = lqProbInf{3};
B = lqProbInf{4};
alpha = lqProbInf{5};
Gamma = lqProbInf{6};
Sigma = lqProbInf{7};

n = size(A, 1);
m = size(B, 2);
S = size(alpha,1);
kappa = size(alpha,2);

%Cost-to-go linear quadratic weightings and constant term
W = zeros(n,n,kappa, S);
W_s = zeros(n,n,kappa,S);
W(:,:,1,:) = repmat(Q, [1 1 S]);

W_s_sum_prev = -1*eye(n);

threshold = 0.000001;


%% SOLUTION
while 1
    for j = 1:kappa
        for s = 1:S
            W_s(:,:,j,s) = W(:,:,j,s);
            % Higher order statistics for uncertainty source s
            if j > 1
                for rc = 1 : j-1
                    W_s(:,:,j,s) = W_s(:,:,j,s) + (j-1) * nchoosek(j-1,rc-1) * W(:,:,rc,s) * Gamma(:,:,s) * Sigma(:,:,s) *Gamma(:,:,s)' * W_s(:,:,j-rc,s);
                end
            end
            
        end
    end
    
    W_s_sum = zeros(n,n);
    % Sum over cumulants
    for j = 1:kappa
        % Average over uncertainty sources
        for s = 1:S
            W_s_sum = W_s_sum + (1/S) * alpha(s,j) * W_s(:,:,j,s);
        end
    end

    % Shortcuts for optimal solution
    G = B'*W_s_sum*A;
    H = R + B'*W_s_sum*B;
    
    [L] = uOptimalFB(G,H);
    
    if  sum(sum(abs(W_s_sum_prev - W_s_sum))) / sum(sum(abs(W_s_sum))) < threshold
        break;
    else
        W_s_sum_prev = W_s_sum;
    end

    % Update cost-to-go
    for j = 1:kappa
        for s = 1:S
            
            if j == 1 % Update Expected cost
                
                % Shortcuts for each marginal computation (each uncertainty source)
                G = B'*W_s(:,:,j,s)*A;
                H = R + B'*W_s(:,:,j,s)*B;
                
                W(:,:,j,s) = Q + A'*W_s(:,:,j,s)*A + ...
                       L'*H*L + L'*G + G'*L;
                
            else % Update Higher order cumulants
                
                % Shortcuts for each marginal computation (each uncertainty source)
                G = B'*W_s(:,:,j,s)*A;
                H = B'*W_s(:,:,j,s)*B;
                
                W(:,:,j,s) =  A'*W_s(:,:,j,s)*A + ...
                       L'*H*L + L'*G + G'*L ;
            end
            
        end 
    end

end

end


function [L] = uOptimalFB(G, H)
    % eigenvalue decomposition, modify eigenvalues
    epsilon = 0;
    [V,D] = eig(H);
    d = diag(D);
    d(d<epsilon) = epsilon;

    % inverse modified Hessian, unconstrained control law
    H1 = V*diag(1./d)*V';
    L = -H1*G;
end

