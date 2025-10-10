%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      cqo1.m
%
%  Purpose:   To demonstrate how to solve a small conic quadratic
%              optimization problem using the MOSEK Matlab Toolbox.
%

function cqo1()

clear prob;

[r, res] = mosekopt('symbcon');
% Specify the non-conic part of the problem.

prob.c   = [0 0 0 1 1 1];
prob.a   = sparse([1 1 2 0 0 0]);
prob.blc = 1;
prob.buc = 1;
prob.blx = [0 0 0 -inf -inf -inf];
prob.bux = inf*ones(6,1);

% Specify the cones as affine conic constraints.
% Two conic constrains: one with QUAD, one with RQUAD, both of dimension 3

prob.accs = [res.symbcon.MSK_DOMAIN_QUADRATIC_CONE 3 res.symbcon.MSK_DOMAIN_RQUADRATIC_CONE 3];

% The matrix such that f * x = [x(4), x(1), x(2), x(5), x(6), x(3)]

prob.f = sparse( 1:6, [4, 1, 2, 5, 6, 3], ones(1, 6) );

% That implies:
%  (x(4), x(1), x(2)) \in QUAD_3
%  (x(5), x(6), x(3)) \in RQUAD_3

% Optimize the problem. 

[r,res]=mosekopt('minimize',prob);

% Display the primal solution.

res.sol.itr.xx'
