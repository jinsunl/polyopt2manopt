% clear all,close all,clc;
% prog = spotsosprog;
% xdim = 1;
% x = msspoly('x',xdim);
% prog = prog.withIndeterminate(x);
% obj = x^2+4*x+5;
% degree = deg(obj,x);
% [prog, lambda] = prog.newFree(1);
% K = [x+3;1-x];
% [prog,tk,f,~] = sosOnK_return_f(prog,obj+lambda,x,K,degree);
% options = spot_sdp_default_options();
% options.verbose = 2;
% sol = prog.minimize(lambda, @spot_mosek_sos, options);
% obj_sos = -double(sol.eval(lambda))

%%
clear all,close all,clc;
x = msspoly('x');
lambda = msspoly('l');
f = 3+2*x+x^2;
f = 6+4*x+6*x^2+4*x^3+x^4;
% f = 35-4*x-7*x^2-16*x^3+5*x^4-x^5+x^6;
% f = 3.5+6.3*x-8.7*x^2-5*x^3+6*x^4;
f = 3.5+6.3*x-8.7*x^3-5*x^5+6*x^8;


f_sos = f-lambda;


% x = msspoly('x',3);
% f = x'*[3 -2 1;-2 5 -3;1 -3 6]*x+6;
% f_sos = f-lambda;

x = msspoly('x',10);
vec = monomials(x,0:2);
temp = randn(length(vec),round(length(vec)/2));
temp = temp*temp';
f = vec'*temp*vec+8;
f_sos = f-lambda;


h = [x-5];
% h = [];
if isempty(h),
        dh = 0;
    else
        dh = full(deg(h,x));
end
i = round(max(full(deg(f,x)),dh)/2);

[ F,A,C,Mandim ,B,b] = pre_process( f,h,x,i );

% %%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(A),
%     A{i} = 0*A{i};
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%% parameters %%%%%%%%%%%%%%%
Rank = 1;                   %
N =  10000000;                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

manifold = symfixedrankYYfactory(Mandim,Rank);
problem.M = manifold;

% generate cost function, gradient and hessian
% part1: cost function
problem.cost = @(Y)cost( Y,F,A,C,N,B,b );

% part2: gradient
problem.egrad = @(Y)egrad(Y,F,A,C,N,B,b);

% % part3: Hessian
% problem.ehess = @(Y,U)ehess(Y,U,F,A,C,N,B,b);
% 
% % manopt
checkgradient(problem); pause;
checkhessian(problem); pause;

options.maxiter = 100000;

[Y,cost] = trustregions(problem,[],options);
res = Y
cost = Y'*F*Y



























