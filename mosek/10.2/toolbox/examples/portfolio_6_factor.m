%%
%  File : portfolio_6_factor.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 of factor type.
% 
%%
function portfolio_6_factor()

n      = 8;
w      = 1.0;
mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
gammas = [0.24 0.28 0.32 0.36 0.4 0.44 0.48];

B    = [ 0.4256, 0.1869; ...
         0.2413, 0.3877; ...
         0.2235, 0.3697; ...
         0.1503, 0.4612; ...
         1.5325, -0.2633; ...
         1.2741, -0.2613; ...
         0.6939, 0.2372; ...
         0.5425, 0.2116 ];

S_F = [0.0620, 0.0577; ...
       0.0577, 0.0908 ];

theta = [0.0720 0.0508 0.0377 0.0394 0.0663 0.0224 0.0417 0.0459];

%% Compute the small factorization required for the model
P  = chol(S_F)';
G_factor  = B * P;
G_factor_T = G_factor';

disp(sprintf('\nBasing Markowitz portfolio optimization using a factor model'));
for gamma=gammas
    [er, x] = FactorMarkowitz(n,mu,G_factor_T,theta,x0,w,gamma);
    disp(sprintf('Expected return: %.4e Std. deviation: %.4e', mu'*x, gamma));
end

%
%    Purpose:
%        Computes the optimal portfolio for a given risk 
%     
%    Input:
%        n: Number of assets
%        mu: An n dimmensional vector of expected returns
%        G_factor_T: The factor (dense) part of the factorized risk
%        theta: specific risk vector
%        x0: Initial holdings 
%        w: Initial cash holding
%        gamma: Maximum risk (=std. dev) accepted
%     
%    Output:
%        Optimal expected return and the optimal portfolio     
%
function [er, x] = FactorMarkowitz(n,mu,G_factor_T,theta,x0,w,gamma)

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

% An affine conic constraint: [gamma, G_factor_T*x, sqrt(theta).'x ] in quadratic cone
prob.f = sparse([ zeros(1,n); G_factor_T; diag(sqrt(theta)) ]);
prob.g = [gamma; zeros(size(G_factor_T,1)+n,1)];
prob.accs = [ res.symbcon.MSK_DOMAIN_QUADRATIC_CONE size(G_factor_T,1)+n+1 ];

% Maximize problem and return the objective value
[rcode,res] = mosekopt('maximize echo(0)', prob, []);

% Check if the interior point solution is an optimal point
% See https://docs.mosek.com/latest/toolbox/accessing-solution.html about handling solution statuses.
if startsWith(res.rcodestr, 'MSK_RES_ERR') | ~strcmp(res.sol.itr.solsta, 'OPTIMAL')
   throw(MException('mosek:solsta', 'Not an optimal solution'));
end

x = res.sol.itr.xx;    
er = mu'*x;
