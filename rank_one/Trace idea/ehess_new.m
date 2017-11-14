function [ fhess ] = ehess_new(Y,U,F,A_stack,stackA,CC,N,dimA)
%This function generates the Euclidian Hessian of the cost function. (user
%design)
%   Input: user indecates
    
    Mdot = U*Y'+Y*U';
    M = Y*Y';
    dhM = 0;
    dhMdot = 0;
if ~isempty(A_stack),
    M_tilde = A_stack*M(:);
    CM = M_tilde'*CC;
    dhM = 2*N*DIAGfun(CM,dimA)*stackA;
    
    Mdot_tilde = A_stack*Mdot(:);
    CMdot = Mdot_tilde'*CC;
    dhMdot = 2*N*DIAGfun(CMdot,dimA)*stackA;
end
    fhess = 2*dhM*U + 2*dhMdot*Y + 2*F*U;
end

