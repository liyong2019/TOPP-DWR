%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      opt_server_sync.m
%
%  Purpose : Demonstrates solving a problem remotely using the OptServer.
%
%            url should be 'http://server:port' or 'https://server:port'
%            cert is the path to a TLS certificate, if using https
%

function opt_server_sync(inputfile, url, cert)
clear prob;
clear param;
clear optserver;

% We read some problem from a file or set it up here
cmd = sprintf('read(%s)', inputfile);
[r,res] = mosekopt(cmd);
prob = res.prob;

% OptServer data
optserver.host = url;
param.MSK_SPAR_REMOTE_TLS_CERT_PATH = cert;

% Perform the optimization with full log output.
[r,res] = mosekopt(sprintf('%s echo(10)', prob.objsense), prob, param, [], optserver); 

% Use the optimal x solution.
xx = res.sol.bas.xx;

