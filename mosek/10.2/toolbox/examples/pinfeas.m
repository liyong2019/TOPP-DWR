%%
%  File : pinfeas.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Purpose: Demonstrates how to fetch a primal infeasibility certificate
%           for a linear problem
%%

function pinfeas()

% In this example we set up a simple problem
[prob] = testProblem();

% Perform the optimization.
[rcode, res] = mosekopt('minimize', prob);

% Check problem status
if strcmp(res.sol.itr.prosta, 'PRIMAL_INFEASIBLE') 
    % Set the tolerance at which we consider a dual value as essential
    eps = 1e-7;

    disp("Variable bounds important for infeasibility: ");
    analyzeCertificate(res.sol.itr.slx, res.sol.itr.sux, eps);
        
    disp("Constraint bounds important for infeasibility: ")
    analyzeCertificate(res.sol.itr.slc, res.sol.itr.suc, eps);
else 
    disp("The problem is not primal infeasible, no certificate to show");
end

% Set up a simple linear problem from the manual for test purposes
function [prob] = testProblem()
prob = [];
prob.c = [1, 2, 5, 2, 1, 2, 1];
prob.a = sparse([1,1,2,2,3,3,3,4,4,5,6,6,7,7],...
                [1,2,3,4,5,6,7,1,5,2,3,6,4,7],...
                [1,1,1,1,1,1,1,1,1,1,1,1,1,1],... 
                7, 7);
prob.blc = [-inf, -inf, -inf, 1100, 200, 500, 500];
prob.buc = [200, 1000, 1000, 1100, 200, 500, 500];
prob.blx = [0, 0, 0, 0, 0, 0, 0];
prob.bux = repmat(inf, 1, 7);
prob

% Analyzes and prints infeasibility contributing elements
% sl - dual values for lower bounds
% su - dual values for upper bounds
% eps - tolerance for when a nunzero dual value is significant
function analyzeCertificate(sl, su, eps)
n = size(sl);
for i=1:n 
    if abs(sl(i)) > eps 
        disp(sprintf("#%d: lower, dual = %e", i, sl(i))); 
    end
    if abs(su(i)) > eps 
        disp(sprintf("#%d: upper, dual = %e", i, su(i)));
    end
end

