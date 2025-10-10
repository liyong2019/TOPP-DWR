%%
%  File : portfolio_1_basic.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :
%    Implements a basic portfolio optimization model.
% 
%%
function portfolio_1_basic()

n      = 8;
w      = 59.0;
mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
x0     = [8.0 5.0 3.0 5.0 2.0 9.0 3.0 6.0]';
gamma  = 36.0;
GT     = [  0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638; ...
            0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506; ...
            0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914; ...
            0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742; ...
            0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ];

disp('Basic Markowitz portfolio optimization')
[er, x] = BasicMarkowitz(n,mu,GT,x0,w,gamma);
disp(sprintf('Expected return: %.4e Std. deviation: %.4e', er, gamma));


%
%    Purpose:
%        Computes the optimal portfolio for a given risk 
%     
%    Input:
%        n: Number of assets
%        mu: An n dimmensional vector of expected returns
%        GT: A matrix with n columns so (GT')*GT  = covariance matrix
%        x0: Initial holdings 
%        w: Initial cash holding
%        gamma: Maximum risk (=std. dev) accepted
%     
%    Output:
%        Optimal expected return and the optimal portfolio     
%
function [er, x] = BasicMarkowitz(n,mu,GT,x0,w,gamma)

[rcode, res] = mosekopt('symbcon echo(0)');
prob = [];

% Objective vector - expected return
prob.c = mu;

% The budget constraint  - e'x = w + sum(x0)
prob.a = ones(1,n);
prob.blc = w + sum(x0);
prob.buc = w + sum(x0);

% Bounds exclude shortselling
prob.blx = zeros(n,1);
prob.bux = inf*ones(n,1); 

% An affine conic constraint: [gamma, GT*x] in quadratic cone
prob.f = sparse([ zeros(1,n); GT ]);
prob.g = [gamma; zeros(n,1)];
prob.accs = [ res.symbcon.MSK_DOMAIN_QUADRATIC_CONE n+1 ];

% Maximize problem and return the objective value
[rcode,res] = mosekopt('maximize echo(0)', prob, []);

% Check if the interior point solution is an optimal point
% See https://docs.mosek.com/latest/toolbox/accessing-solution.html about handling solution statuses.
if startsWith(res.rcodestr, 'MSK_RES_ERR') | ~strcmp(res.sol.itr.solsta, 'OPTIMAL')
   throw(MException('mosek:solsta', 'Not an optimal solution'));
end

x = res.sol.itr.xx;    
er = mu'*x;
