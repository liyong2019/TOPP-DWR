%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      mindisk.m
%
%  Purpose:    A mixed-integer conic quadratic
%              problem.
%
%  Given points (x(i),y(i)) and a bound L it finds 
%  the smallest disk which covers some L of the points.
%%

function mindisk()
% Coordinates of the given points (example)
x = [5.108736e+02 4.848650e+02 7.847473e+01 4.703155e+02 5.461198e+02 5.690343e+02 8.041608e+02 7.777201e+02 4.740465e+02];
y = [-1.906827e+02 -6.262071e+02 -6.383325e+02 -7.925310e+02 -3.453408e+02 -7.685773e+02 -9.449286e+02 -1.694813e+02 -7.486016e+02];

% How many points at least should be covered
L = 5;

% Big-M for the interger formulation
M = 3000;

% Initialize symbolic constants
[r,res] = mosekopt('symbcon');
sym = res.symbcon;

% Variables [R x y u1 ... un]
n = length(x);
prob = [];
prob.c = [1; zeros(n+2,1)];

% Constraint sum(u) >= L
prob.a = sparse([0 0 0 ones(1,n)]);
prob.blc = [L];
prob.buc = [inf];

% ui are binary
prob.ints.sub = 4:n+3;
prob.blx = [-inf;-inf;-inf; zeros(n,1)];
prob.bux = [ inf; inf; inf; ones(n,1)];

% A sequence of constraints 
% (R + M(1-ui), x-xi, y-yi) \in Q3
prob.f=[];
prob.g=[];
for i=1:n
	prob.f=[prob.f; eye(3) sparse(1, i, -M, 3, n)];
	prob.g=[prob.g; M; -x(i); -y(i)];
end
prob.accs = repmat([sym.MSK_DOMAIN_QUADRATIC_CONE 3], 1,n);

% Solve
[r,res] = mosekopt('minimize', prob);
sol = res.sol.int.xx;
R = sol(1);
X = sol(2);
Y = sol(3);

% Optionally plot the result
%scatter(x, y, '.');
%hold on;
%fplot(@(t) X + R*cos(t), @(t) Y + R*sin(t))
%axis equal;
end
