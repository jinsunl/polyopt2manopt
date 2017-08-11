function [ f ] = cost( Y,AA,M,const,N,J )
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates

    X = Y*Y';
    b = J*reshape(X,[],1)-const;
    lambda = AA*b;
    temp = M*lambda-b;
    f = -lambda(1)+N*norm(temp,'fro')^2;


end

