function [ junk ] = Aug_test( x,vec,F,N,stackA,A_stack,dimA)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
     lengthA = size(stackA,1)/dimA;
     y = double(subs(vec,x,randn(length(x),1)));
     M = y*y';
     Mtilde = A_stack*M(:);
     UpperLeft = 2*F+2*N*DIAGfun(Mtilde',dimA)*stackA;
     UpperRight = [];
     for i = 1:lengthA,
         UpperRight = [UpperRight, stackA((i-1)*dimA+1:i*dimA,:)*y];
     end
     AugMat = [UpperLeft, UpperRight;UpperRight', -1/2/N*eye(size(UpperRight,2))];
     junk = eig(AugMat);

end

