function [x, u, L, cost, epsilon] = iterativeLQSolver(fnDyn, fnCost, fnOpt, dt, n, x0, u0, varargin)
%  This implementation builds upon the code provided by E. Todorov. from
%  the paper
%
%  Todorov, E. and Li, W. (2005) A generalized iterative LQG method for
%  locally-optimal feedback control of constrained nonlinear stochastic
%  systems. In proceedings of the American Control Conference, pp 300-306
%
%  Although mostly modified, the main structure of the implementation holds.
%  If you use this code please cite 
%
%  Medina, J.R. and Hirche, S. (2015) Uncertainty-dependent Locally-Optimal 
% Control for Robot Control Considering High-order Cost Statistics
%
%
% [x_, u_, L, cost ] = 
%       iterativeLQSolver(fnDyn, fnCost, dt, n, x0, u0, varargin)
%
% Iterative LQ control for nonlinear plants with bounded controls
%
%DYNAMICS:  x(t+1) = x(t) + dt*fnDyn(x(t), u(t))
%               where x(1) = x0, t = 1:n-1
%COST:      fnCost(x(n)) + sum_t dt*fnCost(x(t),u(t),t)
%CONTROL:   u(t) = u_(t) + L(:,:,t) * (x(t) - x_(t))
%
%INPUTS:
% fnDyn - handle of dynamics function, which must be in the format:
%           [xdot, xdot_x, xdot_u, Gamma, Sigma] = fnDyn(x, u)
%          where x is state, u is control, xdot is time derivative of x,
%          xdot_x and xdot_u are partial derivatives of xdot
%          if only xdot is required, x and u may be matrices
% fnCost- handle of cost function, which must be in the format:
%           [l, l_x, l_xx, l_u, l_uu, l_ux] = fnCost(x, u, t)
%          where l is the cost, and l_x etc. are its partial derivatives
%          t is the current time step; t=NaN means final cost
% fnOpt- handle of optimizer function. In our case either an cost
%        cumulant kcc or a risk-sensitive leqr solver. It takes as an input
%        an lqProb structure.
% dt    - time step for discrete-time simulation
% n     - number of time steps
% x0    - initial state
% u0    - initial control, either vector or matrix with n-1 columns; if u0
%          is vector, it is repeated n-1 times to form a control sequence
% theta - weights for the cumulants
% uMin  - [optional] minimum control, scalar or vector (-Inf: unbounded)
%          if scalar, the bound is the same for all elements of u
% uMax  - [optional] maximum control, scalar or vector (+Inf: unbounded)
%
%OUTPUTS:
% x    - average trajectory under optimal controller (n column matrix)
% u    - open-loop component of optimal controller (n-1 column matrix)
% L     - time-varying feedback gains of optimal controler (3D array)
% cost  - expected cost under optimal controller
% epsilon - lambda value for the line search

global flgFix;

%---------------------- user-adjustable parameters ------------------------

epsilonMin = 1e-10;       % exit if lambda exceeds threshold
relConverge = 1e-5;    % exit if relative improvement below threshold
flgPrint = 1;           % show cost- 0:never, 1:every iter, 2:final
maxValue = 1E+15;       % upper bound on states and costs (Inf: none)
epsilonInit = 1;


%---------------------------------------------- get optional arguments
if nargin>=8
    uMin = varargin{1};
else
    uMin = -Inf;
end
if nargin>=9
    uMax = varargin{2};
else
    uMax = Inf;
end
if nargin>=10
    maxIter = varargin{3};
else
    maxIter = 1000;
end
if nargin>=11
    fnPlot = varargin{4};
    doPlot = true;
    figure;
else
    doPlot = false;
end

szX = size(x0, 1);          % size of state vector
szU = size(u0, 1);          % size of control vector



%% -------------------- INITIALIZATION

L = zeros(szU, szX, n-1);   % init feedback gains

if size(u0,2)==1            % init control sequence
    u = repmat(u0, [1,n-1]);
else
    u = u0;
end

if isscalar(uMin)           % transform scalar arguments into vectors
    uMin = repmat(uMin, [szU,1]);
end
if isscalar(uMax)
    uMax = repmat(uMax, [szU,1]);
end


flgFix = 0;                 % clear large-fix global flag

% In the following lqProb is a cell with :
% 1- Q                     Quadratic State Weighting matrix - n x n x T
% 2- R                     Quadratic Control input weighting - m x m x T-1
% 3- q                     Linear State Weighting matrix - n x T
% 4- r                     Linear Control input weighting - m x T-1
% 5- q0                    Constant cost term - T
% 6- A                     System matrix - n x n x T , where n is the dimension of the state
% 7- B                     Control matrix - n x m x T , where m is the control input dimension
% 8- Gamma                 Control noise matrix - n x qs x S noise inputs
% 9- Sigma                 covariance matrix - qs x qs x T x S
% 10- x0                    Initial state
% (Optional parameters for iterative solution only -> Problem defined on deviations)
% 11- u_nom                 Nominal trajectory m x T-1
% 12- uMin                  Min input m x 1
% 13- uMax                  Max input m x 1


%% ------------------ 2 STEP - OPTIMIZATION LOOP    1st) FF only   2nd) FF + FB
for m = 1:2
    
    %------ STEP 1: Given initial nominal u, simulate dynamics and linearize/quadratize the problem around 
    %-------------  the resulting trajectory to get expected cost, trajectory and local LQ problem
    
    [x, cost, lqProb] = simulate(fnDyn, fnCost, fnOpt, dt, x0, u, L, maxValue);
    
    epsilon = epsilonInit;
    finished = false;
    

    for iter = 1:maxIter 
    %% MAIN ITERATION

        %------ STEP 2: compute optimal control law and cost -> Solve kcc problem
        % Set constraints and nominal trajectory of the lq problem for optimal solution
        lqProb{11} = u;
        lqProb{12} = uMin;
        lqProb{13} = uMax;

        % Solve local LQ problem and get optimal control deviations l and L
        [l, L, ~, ~, ~] = fnOpt(lqProb);


        while 1
        %% LINE SEARCH
            %------ STEP 3: Obtain new control sequence considering control input limits and line search parameter epsilon
            dx = zeros(szX,1);
            unew = zeros(szU, n-1);
            A = lqProb{6};
            B = lqProb{7};
            for k=1:n-1
                du = epsilon*(l(:,k)) + L(:,:,k)*dx;
                du = min(max(du+u(:,k),uMin),uMax) - u(:,k);
                dx = A(:,:,k)*dx + B(:,:,k)*du;  
                unew(:,k) = u(:,k) + du;
            end


            
            %----- Simulate dynamics and linearize/quadratize the problem around resulting trajectory
            %----- to approximate cost (in the first run only FF)
            if m == 1
                % Consider only Feedforward component unew
                [xnew, costtogonew, new_lqProb] = simulate(fnDyn, fnCost, fnOpt, dt, x0, unew, zeros(size(L)), maxValue);
            else
                % Consider also Feedback unew + L \delta x
                [xnew, costtogonew, new_lqProb] = simulate(fnDyn, fnCost, fnOpt, dt, x0, unew, L, maxValue);
            end
            costnew = costtogonew(1);

            
            
            % Plot the current solution
            if (doPlot)
                fnPlot(xnew, unew, L);
            end
            % Print
            if flgPrint==1
                fprintf('Iteration = %d;  Cost = %.4f; MinCost = %.4f; logEps = %.1f\n', ...
                    iter, costnew, cost, log10(epsilon) );
            end


            
            %------ STEP 4: Check for acceptance and update nominal trajectory and LQ problem
            if (costnew< cost) && (~isinf(costnew))

                % Set new trajectory
                u = unew;
                x = xnew;
                
                % Set the previous Linearization/quadratization as the current one
                lqProb = new_lqProb;
                

                % Reset epsilon for the next iteration
                epsilon = epsilonInit;

                
                % Check for convergence
                if iter>1 && ((abs(costnew - cost)/abs(costnew)) < relConverge)
                    display('improvement too small. Done...')
                    cost = costnew;
                    
                    % FINAL CONDITION -> Stop line search and main iteration
                    finished = true;
                    break;          
                    
                end
                cost = costnew;
                
                % Cost decreased -> stop line search
                break;
                
            else
                % Decrease epsilon for line search
                epsilon = epsilon / 2;


                % Check for 'convergence'
                if epsilon<epsilonMin
                    display('epsilon<epsilonMin')
                    
                    % FINAL CONDITION -> Stop line search and main iteration
                    finished = true;
                    break;
                end
            end
        end

        if finished
            break;
        end
    end
end

if flgFix
    warning('ikcc had to limit large numbers, results may be inaccurate');
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulate controlled system, compute trajectory and cost
function [x, cost, lqProb] = simulate(fnDyn, fnCost, fnOpt, dt, x0, u, L, maxValue)
    % get sizes
    szX = size(x0,1);
    szU = size(u,1);
    n = size(u,2)+1;

    % initialize simulation
    x = zeros(szX, n);
    x(:,1) = x0;
    
    Q = zeros(szX,szX,n);
    R = zeros(szU,szU,n-1);
    q = zeros(szX,n);
    r = zeros(szU,n-1);
    q0 = zeros(n);
    A = zeros(szX,szX,n);
    B = zeros(szX,szU,n);

    % run simulation with substeps and approximate dynamics and cost 
    % along new trajectory
    for k = 1:n-1
        [xdot, f_x, f_u, Gamma_k, Sigma_k] = fnDyn(x(:,k), u(:,k));
        x(:,k+1) = fixBig(x(:,k) + dt * xdot, maxValue);
        A(:,:,k) = eye(szX) + dt * f_x;
        B(:,:,k) = dt * f_u;
        Gamma(:,:,k,:) = sqrt(dt)*Gamma_k;
        Sigma(:,:,k,:) = Sigma_k;

        [l0,l_x,l_xx,l_u,l_uu] = fnCost(x(:,k), u(:,k), k);
        q0(k) = dt * l0;
        q(:,k) = dt * l_x;
        Q(:,:,k) = dt * l_xx;
        r(:,k) = dt * l_u;
        R(:,:,k) = dt * l_uu;
    end
    [q0(n),q(:,n),Q(:,:,n)] = fnCost(x(:,n), zeros(szU,1), NaN);  % final cost
    
    lqProb = {Q, R, q, r, q0, A, B, Gamma, Sigma, zeros(szX,1)};

    % Compute expected cost
    [~, ~, costtogo, ~, ~] = fnOpt(lqProb, zeros(size(u)), L);
    cost = costtogo(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Limit numbers that have become suspiciously large
function s = fixBig(s, maxValue)
    global flgFix;

    ibad = (abs(s)>maxValue);
    s(ibad) = sign(s(ibad)) * maxValue;

    if any(ibad(:))
        flgFix = 1;
    end
end