%%
%  File : portfolio_2_frontier.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model.
%                 Determines points on the efficient frontier.
% 
%%

function portfolio_2_frontier()

% Change this to produce the plots
usegraphics = false;

n      = 8;
w      = 1.0;
mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
gamma  = 36.0;
GT     = [  0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638; ...
            0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506; ...
            0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914; ...
            0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742; ...
            0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ];

% Some predefined alphas are chosen
alphas = [0.0 0.01 0.1 0.25 0.30 0.35 0.4 0.45 0.5 0.75 1.0 1.5 2.0 3.0 10.0];

front = EfficientFrontier(n,mu,GT,x0,w,alphas);
disp(sprintf('\nEfficient frontier'));
disp(sprintf('%-12s  %-12s  %-12s', 'alpha', 'return', 'risk'));
for p = front.'
    disp(sprintf('%-12.4f  %-12.4e  %-12.4e', p(1), p(2), p(3)));
end

if usegraphics
    plot(front(:,3),front(:,2));
end


%
%    Purpose:
%        Computes several portfolios on the optimal portfolios by 
%
%            for alpha in alphas: 
%                maximize   expected return - alpha * variance
%                subject to the constraints  
%        
%    Input:
%        n: Number of assets
%        mu: An n dimmensional vector of expected returns
%        GT: A matrix with n columns so (GT')*GT  = covariance matrix
%        x0: Initial holdings 
%        w: Initial cash holding
%        alphas: List of the alphas
%                    
%    Output:
%        The efficient frontier as list of tuples (alpha, expected return, variance)
%
function frontier = EfficientFrontier(n,mu,GT,x0,w,alphas)

frontier = [];
[rcode, res] = mosekopt('symbcon echo(0)');
prob = [];

% The budget constraint in terms of variables [x; s]
prob.a = [ones(1,n), 0.0];
prob.blc = w + sum(x0);
prob.buc = w + sum(x0);

% No shortselling
prob.blx = [zeros(n,1); -inf];
prob.bux = inf*ones(n+1,1); 

% An affine conic constraint: [s, 0.5, GT*x] in rotated quadratic cone
% In matrix form
% [ 0  1]  [ x ]     [ 0   ]
% [ 0  0]  [   ]  +  [ 0.5 ]  \in Q_r
% [ GT 0]  [ s ]     [ 0   ]  
prob.f = sparse([ [zeros(1,n), 1.0]; zeros(1, n+1); [GT, zeros(n,1)] ]);
prob.g = [ 0; 0.5; zeros(n, 1) ]
prob.accs = [ res.symbcon.MSK_DOMAIN_RQUADRATIC_CONE n+2 ];

for alpha = alphas
    % Objective mu'*x - alpha*s
    prob.c = [mu; -alpha];

    [rcode,res] = mosekopt('maximize echo(0)',prob,[]);

    % Check if the interior point solution is an optimal point
    % See https://docs.mosek.com/latest/toolbox/accessing-solution.html about handling solution statuses.
    if startsWith(res.rcodestr, 'MSK_RES_ERR') | ~strcmp(res.sol.itr.solsta, 'OPTIMAL')
       throw(MException('mosek:solsta', 'Not an optimal solution'));
    end

    x = res.sol.itr.xx(1:n);
    s = res.sol.itr.xx(n+1);    
       
    frontier = [frontier; [alpha, mu'*x, sqrt(s)] ];
end
