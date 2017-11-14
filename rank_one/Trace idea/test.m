clear all,close all,clc;

% lambda = msspoly('l');
% 
% x = msspoly('x');
% f = 3+2*x+x^2;
% f = (x+0.5)^4+3;
% f = 35-4*x-7*x^2-16*x^3+5*x^4-x^5+x^6;
% f = 3.5+6.3*x-8.7*x^2-5*x^3+6*x^4;
% f = 3.5+6.3*x-8.7*x^3-5*x^5+6*x^8;
% f_sos = f-lambda;


x = msspoly('x',3);
vec = monomials(x,0:2);
temp = randn(length(vec),round(length(vec)/2));
temp = temp*temp';
% load temp
f = vec'*temp*vec+8;

h = [1-x;x+1];
% h = [];





%% spotless
prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,lambda] = prog.newFree(1);
prog = prog.withSOS(f-lambda);
options = spot_sdp_default_options();
options.verbose = 2;
sol = prog.minimize(-lambda, @spot_mosek_sos, options);
res_spot = sol.eval(lambda)


%% manifold
if isempty(h),
        dh = 0;
    else
        dh = full(deg(h,x));
end
i = round(max(full(deg(f,x)),dh)/2);

[ F,A,C,Mandim ,~] = pre_process_new( f,h,x,i,dim );

% more pre-process needed at this mome
temp_idx = find(C{end}==1);
for i = 1:length(C),
    if length(C{i})>temp_idx+1,
        C{i}(temp_idx+2:end) = [];
    end
    A{i} = 0.5*(A{i}+A{i}');
end
% C{end+1} = 0*C{1};
% C{end}(end) = 1;
% A{end+1} = zeros(Mandim);
% A{end}(1,1) = 1;

F = 0.5*(F+F');
% b = C{end};
b = 0*C{1};

F = temp;F(1,1)=F(1,1)+8;

%%%%%%%%% parameters %%%%%%%%%%%%%%%
Rank = 1;                          %
N =  100000000;                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

manifold = euclideanfactory(Mandim,Rank);
manifold = modefied_euclideanfactory(Mandim,Rank);
% manifold = symfixedrankYYfactory(Mandim,Rank);
problem.M = manifold;

% generate cost function, gradient and hessian
% part1: cost function
problem.cost = @(Y)cost( Y,F,A,C,N,b );

% part2: gradient
problem.egrad = @(Y)egrad(Y,F,A,C,N,b);

% % part3: Hessian
problem.ehess = @(Y,U)ehess(Y,U,F,A,C,N,b);
% 
% % manopt
checkgradient(problem); pause;
checkhessian(problem); pause;

options.maxiter = 100000;

[Y,cost] = trustregions(problem,[],options);
res = sign(Y(1))*Y
% cost = Y'*F*Y
% res_spot




























