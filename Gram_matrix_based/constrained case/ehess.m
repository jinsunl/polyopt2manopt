function [ f_hess ] = ehess(Y,U,M,AA,N,B1,B2,I,e,P1,P2,J)
%This function generates the Euclidian Hessian of the cost function. (user
%design)
%   Input: user indecates
    n = length(B1);
    
    f_hess = 0;
    for i = 1:n,
        for j = 1:n,
            A = B2{j}'*J'*(AA'*M'-I)*(M*AA-I)*J*B2{i};
            B = B1{i}*e*e'*B1{j}';
            f_hess = f_hess + A*U*Y'*B*Y+B*Y*U'*A*Y+A*Y*Y'*B*U;
            f_hess = f_hess + B'*U*Y'*A'*Y+A'*Y*U'*B'*Y+B'*Y*Y'*A'*U;
        end
    end
    f_hess = 2*N*f_hess - P1*U + P2*U;

end

