function [l, L, cost, x, u] = leqr_multi(theta, lqProb, varargin)


% -- INPUTS
% LqProb is a cell with:
% - Q                     Quadratic State Weighting matrix - n x n x T
% - R                     Quadratic Control input weighting - m x m x T-1
% - q                     Linear State Weighting matrix - n x T
% - r                     Linear Control input weighting - m x T-1
% - q0                    Constant cost term - T
% - A                     System matrix - n x n x T , where n is the dimension of the state
% - B                     Control matrix - n x m x T , where m is the control input dimension
% - theta                 Risk-sensitivities S
% - Gamma                 Control noise matrix matrix elements - n x qs x S noise inputs
% - Sigma                 covariance matrix - qs x qs x T x S
% - x0                    Initial state
% - (Optional parameters for iterative solution -> Problem defined on state and control deviations)
% - u_nom                 Nominal trajectory m x T-1
% - uMin                  Min input m x 1
% - uMax                  Max input m x 1

% Only for computing cost
% l                     feedforward control component - m x T-1
% L                     feedback control component  - m x n x T-1

% theta                 risk sensitivities S

% -- OUTPUTS
% l                     feedforward control component - m x T-1
% L                     feedback control component  - m x n x T-1
% cost                cost-to-go of the LEQ problem for initial state x0 where u = l + L * x
% x                   simulated trajectory  - n x T
% u                   simulated control input u = l + L * x

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
S = size(theta,1);

% For iterative solutions
if length(lqProb) > 10
    u_nom = lqProb{11};
    uMin = lqProb{12};
    uMax = lqProb{13};
else
    u_nom = zeros(m,T-1);
    uMin = -inf;
    uMax = inf;
end

% If a control law in the form u=l+L*x is given then only the cost is
% computed
if size(varargin,2)>1,
    compute_opt_control = false;
    l = varargin{1};
    L = varargin{2};
else
    compute_opt_control = true;
    L = zeros(m,n,T-1);
    l = zeros(m,T-1);
end

% Flag indicating any NaN or Inf
overflow = false;

%Cost-to-go linear quadratic weightings and constant term
W = zeros(n,n,T,S);
w = zeros(n,T,S);
w0 = zeros(T,S);
F = zeros(T,S);

W_s = zeros(n,n,S);
w_s = zeros(n,S);
w0_s = zeros(S,1);

W(:,:,T,:) = repmat(Q(:,:,T), [1 1 S]);
w(:,T,:) = repmat(q(:,T), [1 S]);
w0(T,:) = ones(1,S).*q0(T);
F(T,:) = ones(1,S);

%% SOLUTION
for k = T-1:-1:1
    % LEQR solution for each uncertainty source
    for s = 1:S    
        % Calculate auxiliary variables W_s, w_s and w0_s
        W_s(:,:,s) = W(:,:,k+1,s);
        w_s(:,s) = w(:,k+1,s);
        w0_s(s) = w0(k+1,s);

        invTerm = inv(inv(Sigma(:,:,k,s)) - theta(s) * Gamma(:,:,k,s)' * W(:,:,k+1,s) * Gamma(:,:,k,s));
        determ = det(eye(size(Sigma,1)) - theta(s) * Gamma(:,:,k,s)' * W(:,:,k+1,s) * Gamma(:,:,k,s) * Sigma(:,:,k,s));
        W_s(:,:,s) = W_s(:,:,s) + theta(s) * W(:,:,k+1,s) * Gamma(:,:,k,s) * invTerm *Gamma(:,:,k,s)' * W(:,:,k+1,s);
        w_s(:,s) = w_s(:,s) + theta(s) .* W(:,:,k+1,s) * Gamma(:,:,k,s) * invTerm *Gamma(:,:,k,s)' * w(:,k+1,s);
        w0_s(s) = w0_s(s) + theta(s) * w(:,k+1,s)' * Gamma(:,:,k,s) * invTerm *Gamma(:,:,k,s)' * w(:,k+1,s) - theta(s)^(-1)*0.5 * log(determ);
    end
    
    W_s_av = (1/S).*sum(W_s,3);
    w_s_av = (1/S).*sum(w_s,2);

    % Shortcuts
    g = r(:,k) + B(:,:,k)'*w_s_av;
    G = B(:,:,k)'*W_s_av*A(:,:,k);
    H = R(:,:,k) + B(:,:,k)'*W_s_av*B(:,:,k);
    
    if any(isinf(H(:))) > 0 || any(isinf(G(:))) || any(isinf(g)) ...
           || any(isnan(H(:))) > 0 || any(isnan(G(:))) || any(isnan(g))
       overflow = true;
       break;
    end

    if compute_opt_control
        % Compute optimal solution
        [l(:,k), L(:,:,k)] = uOptimal(g,G,H,u_nom(:,k),uMin,uMax);
    end
        
    
    for s = 1:S
        % Marginal shortcuts
        g = r(:,k) + B(:,:,k)'*w_s(:,s);
        G = B(:,:,k)'* W_s(:,:,s)*A(:,:,k);
        H = R(:,:,k) + B(:,:,k)'*W_s(:,:,s)*B(:,:,k);
        
        % Update cost-to-go
        W(:,:,k,s) = Q(:,:,k) + A(:,:,k)'*W_s(:,:,s)*A(:,:,k) + ...
               L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k);
        w(:,k,s) = q(:,k) + A(:,:,k)'*w_s(:,s) + ...
               L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
        w0(k,s) = q0(k) + w0_s(s) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g;
    end
end

x = zeros(n,T);
u = zeros(m,T);
cost = zeros(T,1);
x(:,1) = x0;

for k=1:T-1
    u(:,k) = l(:,k) + L(:,:,k)*x(:,k);
    x(:,k+1) = A(:,:,k)*x(:,k) + B(:,:,k)*u(:,k);
    cost(k) = 0;
    for s=1:S
        cost(k) = cost(k) + (1/2).*x(:,k)'*W(:,:,k,s)*x(:,k) + x(:,k)'*w(:,k,s) + w0(k,s);
    end
    cost(k) = (1/S).*cost(k);
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
