%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      pow1.m
%
%  Purpose:   To demonstrate how to solve the problem
%
%     maximize x^0.2*y^0.8 + z^0.4 - x
%           st x + y + 0.5z = 2
%              x,y,z >= 0

function pow1()

clear prob;

[r, res] = mosekopt('symbcon');

% Specify the non-conic part of the problem.

% Variables number 1,2,3 correspond to x,y,z, variables 4,5 are auxiliary
prob.c   = [-1 0 0 1 1];
prob.a   = [1 1 0.5 0 0];
prob.blc = [2.0];
prob.buc = [2.0];
prob.blx = [-inf -inf -inf -inf -inf];
prob.bux = [ inf  inf  inf  inf  inf];

% Specify the cones as affine conic constraints.
% Two conic constrains with the power cone, both of dimension 3

prob.accs   = [res.symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE 3 2 0.2 0.8 res.symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE 3 2 0.4 0.6];

% The matrices such that f * x + g = [x(1), x(2), x(4), x(3), 1, x(5)] 

prob.f = sparse( [1, 2, 3, 4, 6], [1, 2, 4, 3, 5], ones(1, 5) );
prob.g = [0 0 0 0 1 0];

% That means
%
% (x(1), x(2), x(4)) \ in PPOW_3(0.2, 0.8)
% (x(3), 1,    x(5)) \ in PPOW_3(0.4, 0.6)
%
% which is equivalent to
%
% |x(4)| <= x(1)^0.2 * x(2)^0.8
% |x(5)| <= x(3)^0.4 

% Optimize the problem. 

[r,res]=mosekopt('maximize',prob);

% Display the primal solution.

res.sol.itr.xx'
