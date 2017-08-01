function [l, L, cost, x, u, cumulants] = kcc_multi(gammas, lqProb, varargin)

% -- INPUTS
% LqProb is a cell with :
% - Q                     Quadratic State Weighting matrix - n x n x T
% - R                     Quadratic Control input weighting - m x m x T-1
% - q                     Linear State Weighting matrix - n x T
% - r                     Linear Control input weighting - m x T-1
% - q0                    Constant cost term - T
% - A                     System matrix - n x n x T , where n is the dimension of the state
% - B                     Control matrix - n x m x T , where m is the control input dimension
% - Gamma                 Control noise matrix matrix elements - n x qs x S noise inputs
% - Sigma                 covariance matrix - qs x qs x T x S
% - x0                    Initial state
% - (Optional parameters for iterative solution -> Problem defined on deviations)
% - u_nom                 Nominal trajectory m x T-1
% - uMin                  Min input m x 1
% - uMax                  Max input m x 1

% Only for computing cost
% - l                     feedforward control component - m x T-1
% - L                     feedback control component  - m x n x T-1

% -- OUTPUTS
% l                     Optimal feedforward control component - m x T-1
% L                     Optimal feedback control component  - m x n x T-1
% cost                  cost-to-go of the LEQ problem for initial state x0 where u = l + L * x
% x                     simulated trajectory  - n x T
% u                     simulated control input u = l + L * x


%% INITIALIZATION
Q = lqProb{1};
R = lqProb{2};
q = lqProb{3};
r = lqProb{4};
q0 = lqProb{5};
A = lqProb{6};
B = lqProb{7};
Gamma = lqProb{8};
Sigma = lqProb{9};
x0 = lqProb{10};

T = size(Q,3);
n = size(A, 1);
m = size(B, 2);
S = size(gammas,1);
kappa = size(gammas,2);

% For iterative constrained solutions
if length(lqProb) > 10
    u_nom = lqProb{11};
    uMin = lqProb{12};
    uMax = lqProb{13};
else
    u_nom = zeros(m,T-1);
    uMin = -inf;
    uMax = inf;
end

compute_opt_control = true;
if size(varargin,2)>1,
    compute_opt_control = false;
    l = varargin{1};
    L = varargin{2};
else
    L = zeros(m,n,T-1);
    l = zeros(m,T-1);
end

% Flag indicating any NaN or Inf
overflow = false;

%Cost-to-go linear quadratic weightings and constant term
W = zeros(n,n,kappa,T, S);
w = zeros(n,kappa,T, S);
w0 = zeros(kappa,T, S);
H_s = zeros(n,n,kappa,T-1, S);

W_s = zeros(n,n,kappa,T,S);
w_s = zeros(n,kappa,T,S);
w0_s = zeros(kappa,T,S);

W(:,:,1,T,:) = repmat(Q(:,:,T), [1 1 S]);
w(:,1,T,:) = repmat(q(:,T), [1 S]);
w0(1,T,:) = repmat(q0(T), [S 1]);

%% SOLUTION
for k = T-1:-1:1
    % Calculate W_s, w_s and w0_s
    for j = 1:kappa
        for s = 1:S
            
            W_s(:,:,j,k+1,s) = W(:,:,j,k+1,s);
            w_s(:,j,k+1,s) = w(:,j,k+1,s);
            w0_s(j,k+1,s) = w0(j,k+1,s);
            H_s(:,:,j,k,s) = (1/2) * j * W(:,:,j,k+1,s) * Gamma(:,:,k,s) * Sigma(:,:,k,s) * Gamma(:,:,k,s)';

            % Higher order statistics at time k for uncertainty source s
            if j > 1
                for rc = 1 : j-1
                    H_s(:,:,j,k,s) = H_s(:,:,j,k,s) + (j-1) * nchoosek(j-1,rc-1) * W(:,:,rc,k+1,s) * Gamma(:,:,k,s) * Sigma(:,:,k,s) * Gamma(:,:,k,s)' * H_s(:,:,j-rc,k,s);
                    W_s(:,:,j,k+1,s) = W_s(:,:,j,k+1,s) + (j-1) * nchoosek(j-1,rc-1) * W(:,:,rc,k+1,s) * Gamma(:,:,k,s) * Sigma(:,:,k,s) *Gamma(:,:,k,s)' * W_s(:,:,j-rc,k+1,s);
                    w_s(:,j,k+1,s) = w_s(:,j,k+1,s) + (j-1) * nchoosek(j-1,rc-1) * W(:,:,rc,k+1,s) * Gamma(:,:,k,s) * Sigma(:,:,k,s) *Gamma(:,:,k,s)' * w_s(:,j-rc,k+1,s);
                    w0_s(j,k+1,s) = w0_s(j,k+1,s) + (j-1) * nchoosek(j-1,rc-1) * w(:,rc,k+1,s)' * Gamma(:,:,k,s) * Sigma(:,:,k,s) *Gamma(:,:,k,s)' * w_s(:,j-rc,k+1,s);
                end
            end
            
        end
    end
    
    W_s_sum = zeros(n,n);
    w_s_sum = zeros(n,1);
    % Sum over cumulants
    for j = 1:kappa
        % Average over uncertainty sources
        for s = 1:S
            W_s_sum = W_s_sum + (1/S) * gammas(s,j) * W_s(:,:,j,k+1,s);
            w_s_sum = w_s_sum + (1/S) * gammas(s,j) * w_s(:,j,k+1,s);
        end
    end

    % Shortcuts for optimal solution
    g = r(:,k) + B(:,:,k)'*w_s_sum;
    G = B(:,:,k)'*W_s_sum*A(:,:,k);
    H = R(:,:,k) + B(:,:,k)'*W_s_sum*B(:,:,k);

    if any(isinf(H(:))) > 0 || any(isinf(G(:))) || any(isinf(g)) ...
           || any(isnan(H(:))) > 0 || any(isnan(G(:))) || any(isnan(g))
       overflow = true;
       break;
    end
    
    if compute_opt_control
        % Compute optimal FF + FB solution
        [l(:,k), L(:,:,k)] = uOptimal(g,G,H,u_nom(:,k),uMin,uMax);
    end

    % Update cost-to-go
    for j = 1:kappa
        for s = 1:S
            
            if j == 1 % Update Expected cost
                
                % Shortcuts for each marginal computation (each uncertainty source)
                g = r(:,k) + B(:,:,k)'*w_s(:,j,k+1,s);
                G = B(:,:,k)'*W_s(:,:,j,k+1,s)*A(:,:,k);
                H = R(:,:,k) + B(:,:,k)'*W_s(:,:,j,k+1,s)*B(:,:,k);
                
                W(:,:,j,k,s) = Q(:,:,k) + A(:,:,k)'*W_s(:,:,j,k+1,s)*A(:,:,k) + ...
                       L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k);
                w(:,j,k,s) = q(:,k) + A(:,:,k)'*w_s(:,j,k+1,s) + ...
                       L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
                w0(j,k,s) = q0(k) + w0_s(j,k+1,s) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g;% + trace(H_s(:,:,j,k,s));
                
            else % Update Higher order cumulants
                
                % Shortcuts for each marginal computation (each uncertainty source)
                g = B(:,:,k)'*w_s(:,j,k+1,s);
                G = B(:,:,k)'*W_s(:,:,j,k+1,s)*A(:,:,k);
                H = B(:,:,k)'*W_s(:,:,j,k+1,s)*B(:,:,k);
                
                W(:,:,j,k,s) =  (A(:,:,k)'*W_s(:,:,j,k+1,s)*A(:,:,k) + ...
                       L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k) );
                w(:,j,k,s) =  A(:,:,k)'*w_s(:,j,k+1,s) + ...
                       L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
                w0(j,k,s) = w0_s(j,k+1,s) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g + (1/j)*trace(H_s(:,:,j,k,s));
            end
            
        end 
    end
end

x = zeros(n,T);
u = zeros(m,T);
cost = zeros(T);
cumulants = zeros(kappa,S,T);
x(:,1) = x0;

% Compute optimal trajectory, control and cost
for k=1:T-1
    u(:,k) = l(:,k) + L(:,:,k)*x(:,k);
    x(:,k+1) = A(:,:,k)*x(:,k) + B(:,:,k)*u(:,k);
    cost(k) = 0;
    
    % Compute cost cumulants and their weighted(theta(kappa,S)) sum
    for j = 1:kappa
        for s = 1:S
            cumulants(j,s,k) = j*((1/2)*x(:,k)'*W(:,:,j,k,s)*x(:,k) + x(:,k)'*w(:,j,k,s) + w0(j,k,s));
            cost(k) = cost(k) + (1/S) * gammas(s,j) * cumulants(j,s,k);
        end
    end
end

if overflow
    cost = Inf*ones(T,1);
end

end


% From Li,Todorov 2005
function [l,L] = uOptimal(g, G, H, u, uMin, uMax)
    % eigenvalue decomposition, modify eigenvalues
    epsilon = 1e-10;
    [V,D] = eig(H);
    d = diag(D);
    d(d<epsilon) = epsilon;

    % inverse modified Hessian, unconstrained control law
    H1 = V*diag(1./d)*V';
    l = -H1*g;
    L = -H1*G;

    % enforce constraints
    l = min(max(l+u,uMin),uMax) - u;

    % modify L to reflect active constraints
    L((l+u==uMin)|(l+u==uMax),:) = 0;
end

