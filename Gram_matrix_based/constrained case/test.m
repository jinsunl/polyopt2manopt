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
f = 35-4*x-7*x^2-16*x^3+5*x^4;
f = 3.5+6.3*x-8.7*x^2-5*x^3+6*x^4;
f = 3.5+6.3*x-8.7*x^3-5*x^5+6*x^8;
h = [x+5;-x];
f_sos = f-lambda;


% x = msspoly('x',3);
% f = x'*[3 -2 1;-2 5 -3;1 -3 6]*x+6;
% f_sos = f-lambda;

% x = msspoly('x',8);
% vec = monomials(x,0:2);
% temp = randn(length(vec),round(length(vec)/2));
% temp = temp*temp';
% f = vec'*temp*vec+8;
% f_sos = f-lambda;



% [ M_tilde,M_const,b_aux,dim ] = poly2mat( f_sos, x, lambda);
[ Mi,Consti,b_auxi, Mandim ] = poly2mat_constr( f_sos,h, x, lambda);

%%%%%%%%% parameters %%%%%%%%%%%%%%%
Rank = Mandim-2;                   %
N =  10000000;                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

manifold = symfixedrankYYfactory(Mandim,Rank);
problem.M = manifold;

% generate cost function, gradient and hessian
% part1: cost function
temp = [Mi{1};Mi{2}];
AA = [(temp'*temp)\temp', Mi{3}'];
M = [Mi{1};Mi{2};Mi{3}];
b_aux = [b_auxi{1};b_auxi{2};b_auxi{3}];
const = [Consti{1};Consti{2};Consti{3}];
J = b_aux;
problem.cost = @(Y)cost( Y,AA,M,const,N,J );

% part2: gradient
II = eye(Mandim);
e1 = zeros(size(AA,1),1);e1(1) = 1;
e = zeros(Mandim,1);e(1) = 1;
I = eye(size(AA,2));
P1 = 0;
P2 = 0;
for i = 1:Mandim,
   Base1{i} = mat_col_switch(II,1,i); 
   Base2{i} = [zeros(Mandim*(i-1),Mandim);eye(Mandim);zeros(Mandim*(Mandim-i),Mandim)];
   P1 = P1 + Base2{i}'*J'*AA'*e1*e'*Base1{i}';
   P2 = P2 + Base2{i}'*J'*(AA'*M'-I)*(M*AA-I)*const*e'*Base1{i}';
end
P1 = P1 + P1';
P2 = -2*N*(P2+P2');

problem.egrad = @(Y)egrad(Y,M,AA,N,Base1,Base2,I,e,P1,P2,J);

% % part3: Hessian
problem.ehess = @(Y,U)ehess(Y,U,M,AA,N,Base1,Base2,I,e,P1,P2,J);
% 
% % manopt
checkgradient(problem); pause;
checkhessian(problem); pause;
Y = trustregions(problem);
res = AA*(J*reshape(Y*Y',size(Y,1)^2,1)-const);
res(1)
M*res+const-J*reshape(Y*Y',size(Y,1)^2,1)




























