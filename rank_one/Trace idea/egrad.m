function [ fgrad ] = egrad(Y,F,A,C,N,b)
%This function generates the Euclidian gradient of the cost function. (user
%design)
%   Input: user indecates
    
    M = Y*Y';
    fgrad = 0;
    if ~isempty(A),
        n = size(A,3);
        for i = 1:n,
            for j = 1:n,
                fgrad = fgrad + C(:,:,i)'*C(:,:,j)*(A(:,:,i)*trace(A(:,:,j)*M)+A(:,:,j)*trace(A(:,:,i)*M));
            end
            fgrad = fgrad - 2*C(:,:,i)'*b*A(:,:,i);
        end
    end
    fgrad = fgrad*N*2*Y+F*2*Y;

end

