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
% f = [-150 34 65]*[1;x;x^2];
f = 6+4*x+6*x^2+4*x^3+x^4;
f = 35-4*x-7*x^2-16*x^3+5*x^4;
f = 3.5+6.3*x-8.7*x^2-5*x^3+6*x^4;
f = 3.5+6.3*x-8.7*x^3-5*x^5+6*x^8;
f_sos = f-lambda;
% [ M,Basis,Mono ] = fun2quadform( f_sos, x);
% M{1} = M{1}+[0 0 1; 0 -2 0;1 0 0];
% [ M_tilde,M_const ] = quadmat2vec( M,lambda );

% x = msspoly('x',3);
% f = x'*[3 -2 1;-2 5 -3;1 -3 6]*x+6;
% f_sos = f-lambda;

x = msspoly('x',8);
vec = monomials(x,0:2);
temp = randn(length(vec),round(length(vec)/2));
temp = temp*temp';
f = vec'*temp*vec+8;
f_sos = f-lambda;



[ M_tilde,M_const,b_aux,dim ] = poly2mat( f_sos, x, lambda);

% Xdim = round(deg(f)/2)+1;   % wrong for indeterminate variables with multi-demension 
Xdim = dim{1};

Rank = round(Xdim/2);
manifold = symfixedrankYYfactory(Xdim,Rank);
problem.M = manifold;
% generate cost function, gradient and hessian
% part1: cost function
q = 1; % weight in norm matrix when calculating undeterminate solution 
N = 1000000; %regularizer for the constraint of making matrix on the manifold
    [row,col] = size(M_tilde{1});
    if rank(M_tilde{1})>=col,
        AA = (M_tilde{1}'*M_tilde{1})\(M_tilde{1}');
    else
        Q = eye(size(M_tilde{1},2)); 
        Q(1,1) = q;
        Qinv = Q\eye(size(M_tilde{1},2));
        AA = Qinv*M_tilde{1}'*(inv(M_tilde{1}*Qinv*M_tilde{1}'));
    end
    
    
    const = M_const{1};
    J = b_aux{1};
problem.cost = @(Y)cost( Y,AA,M_tilde{1},const,N,J );

% part2: gradient
% generate two basis and other auxiliary parameters
II = eye(Xdim);
e1 = zeros(size(AA,1),1);e1(1) = 1;
e = zeros(Xdim,1);e(1) = 1;
I = eye(size(AA,2));
P1 = 0;
P2 = 0;
for i = 1:Xdim,
   Base1{i} = mat_col_switch(II,1,i); 
   Base2{i} = [zeros(Xdim*(i-1),Xdim);eye(Xdim);zeros(Xdim*(Xdim-i),Xdim)];
   P1 = P1 + Base2{i}'*J'*AA'*e1*e'*Base1{i}';
   P2 = P2 + Base2{i}'*J'*(AA'*M_tilde{1}'-I)*(M_tilde{1}*AA-I)*const*e'*Base1{i}';
end
P1 = P1 + P1';
P2 = -2*N*(P2+P2');

problem.egrad = @(Y)egrad(Y,M_tilde{1},AA,N,Base1,Base2,I,e,P1,P2,J);

% % part3: Hessian
problem.ehess = @(Y,U)ehess(Y,U,M_tilde{1},AA,N,Base1,Base2,I,e,P1,P2,J);
% 
% % manopt
checkgradient(problem); pause;
checkhessian(problem); pause;
Y = trustregions(problem);
res = AA*(J*reshape(Y*Y',size(Y,1)^2,1)-const);
res(1)
M_tilde{1}*res+const-J*reshape(Y*Y',size(Y,1)^2,1)

