function [ fgrad ] = egrad_new(Y,F,A_stack,stackA,CC,N,dimA)
%This function generates the Euclidian gradient of the cost function. (user
%design)
%   Input: user indecates
    
    M = Y*Y';
    dhM = 0;
    if ~isempty(A_stack),
        M_tilde = A_stack*M(:);
        CM = M_tilde'*CC;
        dhM = 2*N*DIAGfun(CM,dimA)*stackA;
    end
    fgrad = dhM*2*Y+F*2*Y;

end

