%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      ceo1.m
%
%  Purpose:   To demonstrate how to solve a small conic exponential
%              optimization problem using the MOSEK Matlab Toolbox.
%

function ceo1()

clear prob;

[r, res] = mosekopt('symbcon');
% Specify the non-conic part of the problem.

prob.c   = [1 1 0];
prob.a   = sparse([1 1 1]);
prob.blc = 1;
prob.buc = 1;
prob.blx = [-inf -inf -inf];
prob.bux = [ inf  inf  inf];

% Specify the affine conic constraint with one exponential cone.

prob.accs   = [res.symbcon.MSK_DOMAIN_PRIMAL_EXP_CONE 3];
prob.f = speye(3);

% prob.accs the domain types, in this case a single exponential cone
% The matrix f is the ientity, meaning that 
%
% I * x \in EXP
%
% which is exactly
%
%   x(1) >= x(2)*exp(x(3)/x(2))

% Optimize the problem. 

[r,res]=mosekopt('minimize',prob);

% Display the primal solution.

res.sol.itr.xx'
