%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      helloworld.m
%
%  The most basic example of how to get started with MOSEK.

prob.a  =  sparse(0,1)  % 0 linear constraints, 1 variable
prob.c  =  [1.0]'       % Only objective coefficient
prob.blx=  [2.0]'       % Lower bound(s) on variable(s)
prob.bux=  [3.0]'       % Upper bound(s) on variable(s)

% Optimize
[r, res] = mosekopt('minimize', prob);

% Print answer
res.sol.itr.xx
