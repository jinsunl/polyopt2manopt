function [ fhess ] = ehess(Y,U,F,A,C,N,b)
%This function generates the Euclidian Hessian of the cost function. (user
%design)
%   Input: user indecates
    
    dhM = 0;
    dhMMdot = 0;
    Mdot = U*Y'+Y*U';
    M = Y*Y';
if ~isempty(A),
    n = size(A,3);
    for i = 1:n,
        for j = 1:n,
            dhM = dhM + C(:,:,i)'*C(:,:,j)*( A(:,:,i)*trace(A(:,:,j)*M)+A(:,:,j)*trace(A(:,:,i)*M) );
            dhMMdot = dhMMdot + C(:,:,i)'*C(:,:,j)*( A(:,:,i)*trace(A(:,:,j)*Mdot)+A(:,:,j)*trace(A(:,:,i)*Mdot) );
        end
        dhM = dhM - 2*C(:,:,i)'*b*A(:,:,i);
    end
end
    fhess = 2*N*(dhM*U + dhMMdot*Y) + 2*F*U;
end

