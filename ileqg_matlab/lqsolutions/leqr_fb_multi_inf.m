function [L] = leqr_fb_multi_inf(lqProbInf, varargin)

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
theta = lqProbInf{5};
Gamma = lqProbInf{6};
Sigma = lqProbInf{7};

n = size(A, 1);
S = size(theta,1);

%Cost-to-go linear quadratic weightings and constant term

W_s = zeros(n,n,S);
W = repmat(Q, [1 1 S]);

%% SOLUTION

W_s_av_prev = -1*eye(n);

threshold = 0.000001;

while 1
    % LEQR solution for each uncertainty source
    for s = 1:S    
        % Calculate auxiliary variables W_s, w_s and w0_s
        invTerm = eye(n) - theta(s) * Gamma(:,:,s)' * W(:,:,s) * Gamma(:,:,s) * Sigma(:,:,s);
        W_s(:,:,s) =  invTerm \ W(:,:,s);
    end
    
    W_s_av = (1/S).*sum(W_s,3);

    % Shortcuts
    G = B'*W_s_av*A;
    H = R + B'*W_s_av*B;
    
    [L] = uOptimalFB(G,H);
    
    if  sum(sum(abs(W_s_av_prev - W_s_av))) / sum(sum(abs(W_s_av))) < threshold
        break;
    else
        W_s_av_prev = W_s_av;
    end
    
    for s = 1:S
        G = B'* W_s(:,:,s)*A;
        H = R + B'*W_s(:,:,s)*B;
        
        % Update cost-to-go
        W(:,:,s) = Q + A'*W_s(:,:,s)*A + ...
               L'*H*L + L'*G + G'*L;
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
