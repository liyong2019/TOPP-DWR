%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      sdo_lmi.m
%
%  Purpose :   To solve a problem with an LMI and an affine conic constrained problem with a PSD term
%
%                 minimize    Tr [1, 0; 0, 1]*X + x(1) + x(2) + 1
%
%                 subject to  Tr [0, 1; 1, 0]*X - x(1) - x(2) >= 0
%                             x(1) [0, 1; 1, 3] + x(2) [3, 1; 1, 0] - [1, 0; 0, 1] >> 0
%                             X >> 0
%%
                
function sdo_lmi()

[r, res] = mosekopt('symbcon');
symbcon = res.symbcon;

% The scalar part, as in linear optimization examples
prob.c = [1.0 1.0];
prob.cfix = 1.0;
prob.a = sparse([], [], [], 0, 2);          % 0 constraints, 2 scalar variables
prob.blc = [];                              % Bounds
prob.buc = [];
prob.blx = [-inf, -inf];
prob.bux = [inf,  inf];

prob.f     =  sparse([1, 1, 2, 3, 3, 4],                     ...
                     [1, 2, 2, 1, 2, 1],                     ...
                     [-1, -1, 3, sqrt(2), sqrt(2), 3],       ...
                     4, 2);
prob.g     = [0, -1, 0, -1];

% Specify the affine conic structure
prob.accs  = [symbcon.MSK_DOMAIN_RPLUS 1 symbcon.MSK_DOMAIN_SVEC_PSD_CONE 3];

% Dimensions of PSD variables
prob.bardim = [2];

% Block triplet format specifying the lower triangular part 
% of the symmetric coefficient matrix 'barc':
prob.barc.subj = [1, 1, 1];
prob.barc.subk = [1, 2, 2];
prob.barc.subl = [1, 1, 2];
prob.barc.val  = [1, 0, 1];

% Block triplet format specifying the lower triangular part 
% of the symmetric coefficient matrix 'barF' for the ACC:
prob.barf.subi = [1, 1];
prob.barf.subj = [1, 1];
prob.barf.subk = [1, 2];
prob.barf.subl = [1, 1];
prob.barf.val  = [0, 1];

% Solve with log output
[r, res] = mosekopt('minimize', prob);

% Print the solution
X = zeros(2);
X([1,2,4]) = res.sol.itr.barx;
X = X + tril(X,-1)';

x = res.sol.itr.xx;

X
x

