% clear all,close all,clc
% manifold = euclideanfactory(2,1);
% manifold = symfixedrankYYfactory(2,1);
% % manifold = stiefelfactory(2,2);
% problem.M = manifold;
% A = rand(2);
% A = A+A';
% F = rand(2);
% F = F+F';
% b = rand(1);
% % problem.cost = @(M)trace(A*M)^2+3;
% % problem.egrad = @(M)2*trace(A*M)*A;
% % 
% % problem.ehess = @(M,U)2*[A(1,1)*reshape(A,4,1)'*reshape(U',4,1), A(2,1)*reshape(A,4,1)'*reshape(U',4,1);...
% %                          A(1,2)*reshape(A,4,1)'*reshape(U',4,1), A(2,2)*reshape(A,4,1)'*reshape(U',4,1)];
% % problem.ehess = @(M,U)2*trace(A*U)*A;
% 
% % problem.cost = @(Y)trace(A*Y*Y')^2 + trace(F*Y*Y');
% % problem.egrad = @(Y)4*trace(A*Y*Y')*A*Y + 2*F*Y;
% % problem.ehess = @(Y,U)4*trace(A*(U*Y'+Y*U'))*A*Y+4*trace(A*Y*Y')*A*U + 2*F*U;
% 
% problem.cost = @(Y)norm(trace(A*Y*Y')-b,'fro')^2 + trace(F*Y*Y');
% problem.egrad = @(Y)( 2*A*trace(A*Y*Y')-2*b*A   )*2*Y  + 2*F*Y;
% problem.ehess = @(Y,U) ( 2*A*trace(A*Y*Y')-2*b*A   )*2*U+( 2*A*trace(A*Y*U'+A*U*Y') )*2*Y + 2*F*U;
% 
% checkgradient(problem);pause
% checkhessian(problem); pause
% MM = trustregions(problem);

%%
clear all,close all,clc;
dim = 5;
manifold = euclideanfactory(dim,1);
manifold = symfixedrankYYfactory(dim,1);
problem.M = manifold;

n = 4;
for i = 1:n,
   A{i} = rand(dim);
   A{i} = A{i}+A{i}';
   C{i} = zeros(n,1);
   C{i}(i) = 1;
end
b = rand(n,1);
F = rand(dim);
F = F+F';

problem.cost = @(Y)cost(Y,A,b,C,F);
problem.egrad = @(Y)egrad(Y,A,b,C,F);


checkgradient(problem);pause
checkhessian(problem); pause
MM = trustregions(problem);

function f = cost(Y,A,b,C,F)
    f = trace(F*Y*Y');
    temp = 0;
    for i = 1:length(A),
        temp = temp + C{i}*trace(A{i}*Y*Y');
    end
    f = f+norm(temp-b,'fro')^2;
end

function fgrad = egrad(Y,A,b,C,F)
    temp = 0;
    for i = 1:length(A),
        for j = 1:length(A),
            temp = temp + C{i}'*C{j}*(A{i}*trace(A{j}*Y*Y') + A{j}*trace(A{i}*Y*Y'));
        end
        temp = temp - 2*C{i}'*b*A{i};
    end
    fgrad = (F+temp)*2*Y;
end

function fhess = ehess(Y,U,A,b,C,F)
    M = Y*Y';
    Mdot = Y*U'+U*Y';
    fhess = 0;
    dMf = F;
    dMfMdot = 0;
    for i = 1:length(A),
        for j = 1:length(A),
            dMf = dMf + C{i}'*C{j}*(A{i}*trace(A{j}*Y*Y') + A{j}*trace(A{i}*Y*Y'));
            dMfMdot = dMfMdot + C{i}'*C{j}*(A{i}*trace(A{j}*Mdot) + A{j}*trace(A{i}*Mdot));
        end
        dMf = dMf - 2*C{i}'*b*A{i};
    end
    fhess = dMf*2*U + dMfMdot*2*Y;
end




% %% 
% manifold = euclideanfactory(3,1);
% manifold = stiefelfactory(3,1);
% A = rand(3);
% A = A+A';
% problem.M = manifold;
% problem.cost = @(M)M'*A*M;
% problem.egrad = @(M)2*A*M;
% problem.ehess = @(M,U)2*A*U;
% checkgradient(problem);
% checkhessian(problem);
% 
% 
% 
% %% test hessian
% n = 3;
% A = sym('a',[n,n],'real');
% X = sym('x',[n,n],'real');
% Z = sym('z',[n,n],'real');
% f = trace(A*X)^2;
% for i = 1:n,
%     for  j = 1:n,
%         e{i}{j}=zeros(n);
%         e{i}{j}(i,j) = 1;
%     end
% end
% hess_real = 0;
% U = Z;
% for i = 1:n,
%     for j = 1:n,
%         for ii = 1:n,
%             for jj = 1:n,
%                 hess_real = hess_real + diff(diff(f,X(i,j)),X(ii,jj))*Z(ii,jj)*e{i}{j};
%             end
%         end
%     end
% end
% hess = 2*[A(1,1)*reshape(A,9,1)'*reshape(U',9,1), A(2,1)*reshape(A,9,1)'*reshape(U',9,1), A(3,1)*reshape(A,9,1)'*reshape(U',9,1);
%           A(1,2)*reshape(A,9,1)'*reshape(U',9,1), A(2,2)*reshape(A,9,1)'*reshape(U',9,1), A(3,2)*reshape(A,9,1)'*reshape(U',9,1);
%           A(1,3)*reshape(A,9,1)'*reshape(U',9,1), A(2,3)*reshape(A,9,1)'*reshape(U',9,1), A(3,3)*reshape(A,9,1)'*reshape(U',9,1)];
% hess-hess_real