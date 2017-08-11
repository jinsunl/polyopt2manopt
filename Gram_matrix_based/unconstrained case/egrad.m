function [ f_grad ] = egrad(Y,M,AA,N,B1,B2,I,e,P1,P2,J)
%This function generates the Euclidian gradient of the cost function. (user
%design)
%   Input: user indecates
    
    n = length(B1);
    
    dlambda_dX = P1*Y;
    
    dfro_dX = 0;
    for i = 1:n,
        for j = 1:n,
            dfro_dX = dfro_dX + B2{j}'*J'*(AA'*M'-I)*(M*AA-I)*J*B2{i}*Y*Y'*B1{i}*e*e'*B1{j}';
        end
    end
    dfro_dX = 2*N*(dfro_dX+dfro_dX')*Y+P2*Y;
    
    f_grad = dfro_dX - dlambda_dX;

end

