function [ f ] = cost( Y,AA,M,const,N,J )
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates

    X = Y*Y';
    b = J*reshape(X,size(X,1)^2,1)-const;
%     [row,col] = size(M);
%     if rank(M)>=col,
%         lambda = (M'*M)\(M'*b);
%     else
%         Q = eye(size(x,1));
%         Q(1,1) = q;
%         Qinv = Q\eye(size(x,1));
%         lambda = Qinv*M'*((M*Qinv*M')\b);
%     end
    lambda = AA*b;
    temp = M*lambda-b;
    f = -lambda(1)+N*norm(temp,'fro')^2;


end

