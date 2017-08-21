clear all, close all,clc;
x = msspoly('x',8);
vec = monomials(x,0:2);
% temp = randn(length(vec),round(length(vec)/2));
% temp = temp*temp';
load temp;
f = vec'*temp*vec+8;

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,lambda] = prog.newFree(1);
prog = prog.withSOS(f-lambda);
% prog = sosOnK(prog,f-lambda,x,x-5,4);
options = spot_sdp_default_options();
options.verbose = 2;
sol = prog.minimize(-lambda, @spot_mosek_sos, options);
sol.eval(lambda)