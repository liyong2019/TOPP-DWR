%%
%  File : portfolio_3_impact.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 with transaction costs of order x^(3/2).
% 
%%

function portfolio_3_impact()

n      = 8;
w      = 1.0;
mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
gamma  = 0.36;
GT     = [  0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638; ...
            0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506; ...
            0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914; ...
            0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742; ...
            0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327; ...
            0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ];

% Somewhat arbirtrary choice of m
m = 1.0e-2*ones(n,1);

[x, t] = MarkowitzWithMarketImpact(n,mu,GT,x0,w,gamma,m);
disp(sprintf('\nMarkowitz portfolio optimization with market impact cost'));
disp(sprintf('Expected return: %.4e Std. deviation: %.4e Market impact cost: %.4e', mu'*x, gamma, m'*t));


%
%        Description:
%            Extends the basic Markowitz model with a market cost term.
%
%        Input:
%            n: Number of assets
%            mu: An n dimmensional vector of expected returns
%            GT: A matrix with n columns so (GT')*GT  = covariance matrix
%            x0: Initial holdings 
%            w: Initial cash holding
%            gamma: Maximum risk (=std. dev) accepted
%            m: It is assumed that  market impact cost for the j'th asset is
%               m_j|x_j-x0_j|^3/2
%
%        Output:
%           Optimal expected return and the optimal portfolio     
%
function [x, t] = MarkowitzWithMarketImpact(n,mu,GT,x0,w,gamma,m)

[rcode, res] = mosekopt('symbcon echo(0)');

% unrolled variable ordered as (x, t)
prob = [];
prob.c = [mu; zeros(n,1)];

In = speye(n);
On = sparse([],[],[],n,n);

% Linear part
% [ e' m' ]  * [ x; t ]  =   w + e'*x0
prob.a   = [ ones(1,n), m' ]; 
prob.blc = [ w + sum(x0) ];
prob.buc = [ w + sum(x0) ];

% No shortselling and no other bounds
prob.blx = [ zeros(n,1); -inf*ones(n,1) ];
prob.bux = inf*ones(2*n,1);

% Affine conic constraints representing [ gamma, GT*x ] in quadratic cone
prob.f = sparse([ zeros(1,2*n); [GT On] ]);
prob.g = [gamma; zeros(n,1)];
prob.accs = [ res.symbcon.MSK_DOMAIN_QUADRATIC_CONE n+1 ];

% Extend the affine conic constraints
% with power cones representing t(i) >= |x(i)-x0(i)|^1.5
fi = [];
fj = [];
g  = [];
fv = repmat([1; 1], n, 1);
for k=1:n
    fi = [fi; 3*k-2; 3*k];
    fj = [fj; n+k; k];
    g  = [g; 0; 1; -x0(k)];
end
prob.f = [prob.f; sparse(fi, fj, fv)];
prob.g = [prob.g; g];
prob.accs = [prob.accs repmat([res.symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE, 3, 2, 2.0, 1.0], 1, n) ];

[rcode,res] = mosekopt('maximize echo(0)',prob,[]);

% Check if the interior point solution is an optimal point
% See https://docs.mosek.com/latest/toolbox/accessing-solution.html about handling solution statuses.
if startsWith(res.rcodestr, 'MSK_RES_ERR') | ~strcmp(res.sol.itr.solsta, 'OPTIMAL')
   throw(MException('mosek:solsta', 'Not an optimal solution'));
end

x = res.sol.itr.xx(1:n);
t = res.sol.itr.xx(n+(1:n));

