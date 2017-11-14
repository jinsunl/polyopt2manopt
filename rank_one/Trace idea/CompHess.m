function [ direc_hess, real_hess,y ] = CompHess( x,vec,hess,stackA,C_stack,F,N,dimA )
% This function try to find out the normal hessian of the cost function. 
% Turns out that direc_hess = real_hess = normal hessian. 

    y = double(subs(vec,x,randn(length(x),1)));
    y = randn(length(y),1);
    M = y*y';
    zero = zeros(length(y),1);
    lengthA = size(C_stack,2);

    real_hess = 2*F;
    for i = 1:lengthA,
        ci = C_stack(:,i);
        Ai = stackA((i-1)*dimA+1:i*dimA,:);
        for j = 1:lengthA,
            cj = C_stack(:,j);
            Aj = stackA((j-1)*dimA+1:j*dimA,:);
            real_hess = real_hess + 4*N*ci'*cj*Aj*trace(Ai*M) + 8*N*ci'*cj*Aj*M*Ai;
        end
    end
    
    direc_hess = real_hess*0;
    for i = 1:length(zero),
        U = zero;
        U(i) = 1;
        direc_hess(:,i) = hess(y,U);
    end
    
    
    eig(direc_hess)
    
    
end

