%%
%  File : portfolio_4_transcost.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 with fixed setup costs and transaction costs
%                 as a mixed-integer problem.
% 
%%

function portfolio_4_transcost()

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

f = 0.01*ones(n,1);
g = 0.001*ones(n,1);

[x, z, y] = MarkowitzWithTransactionsCost(n,mu,GT,x0,w,gamma,f,g);
disp(sprintf('\nMarkowitz portfolio optimization with transactions cost'))
disp(sprintf('Expected return: %.4e Std. deviation: %.4e Transactions cost: %.4e', ...
                 mu'*x,gamma,f'*y+g'*z))


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
%            f: If asset j is traded then a fixed cost f_j must be paid
%            g: If asset j is traded then a cost g_j must be paid for each unit traded
%
%        Output:
%           Optimal expected return and the optimal portfolio     
%
function [x, z, y] = MarkowitzWithTransactionsCost(n,mu,GT,x0,w,gamma,f,g)

[rcode, res] = mosekopt('symbcon echo(0)');

% Upper bound on the traded amount
u = w+sum(x0);

% unrolled variable ordered as (x, z, y)
prob = [];
prob.c = [mu; zeros(2*n,1)];
In = speye(n);
On = sparse([],[],[],n,n);

% Linear constraints
% [ e'  g'  f' ]   [ x ]  =   w + e'*x0
% [ I  -I   0  ] * [ z ]  <=  x0
% [ I   I   0  ]   [ y ]  >=  x0
% [ 0   I  -U  ]          <=  0
prob.a   = [ [ones(1,n), g', f']; In -In On; In In On; On In -u*In ]; 
prob.blc = [ w + sum(x0); -inf*ones(n,1); x0; -inf*ones(n,1) ];
prob.buc = [ w + sum(x0); x0; inf*ones(n,1); zeros(n,1) ];

% No shortselling and the linear bound 0 <= y <= 1 
prob.blx = [ zeros(n,1); -inf*ones(n,1); zeros(n,1) ];
prob.bux = [ inf*ones(2*n,1); ones(n,1) ];

% Affine conic constraints representing [ gamma, GT*x ] in quadratic cone
prob.f = sparse([ zeros(1,3*n); [GT On On];  ]);
prob.g = [gamma; zeros(n,1)];
prob.accs = [ res.symbcon.MSK_DOMAIN_QUADRATIC_CONE n+1 ];

% Demand y to be integer (hence binary)
prob.ints.sub = 2*n+(1:n);

[rcode,res] = mosekopt('maximize echo(0)',prob,[]);

% Check if the integer solution is an optimal point
% See https://docs.mosek.com/latest/toolbox/accessing-solution.html about handling solution statuses.
if startsWith(res.rcodestr, 'MSK_RES_ERR') | ~strcmp(res.sol.int.solsta, 'INTEGER_OPTIMAL')
   throw(MException('mosek:solsta', 'Not an optimal solution'));
end

x = res.sol.int.xx(1:n);
z = res.sol.int.xx(n+(1:n));
y = res.sol.int.xx(2*n+(1:n));
