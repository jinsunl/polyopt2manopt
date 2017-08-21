function [ fgrad ] = egrad(Y,F,A,C,N,B,b)
%This function generates the Euclidian gradient of the cost function. (user
%design)
%   Input: user indecates
    
    fgrad = 0;
if ~isempty(A),
    n = length(A);
    for i = 1:n,
        fgrad = fgrad + 2*N*(B'*C{i}*Y'*A{i}*Y + A{i}'*Y*C{i}'*B*Y + A{i}*Y*Y'*B'*C{i} - A{i}'*Y*C{i}'*b - A{i}*Y*b'*C{i});
%         fgrad = fgrad + A{i}*Y*Y'*B'*C{i} + A{i}'*Y*C{i}'*B*Y + B'*C{i}*Y'*A{i}*Y - A{i}*Y*b'*C{i} - A{i}'*Y*C{i}'*b;
        for j = 1:n,
            fgrad = fgrad + 2*N*(A{i}*Y*Y'*A{j}'*Y*C{j}'*C{i}+...
                                 A{i}'*Y*C{i}'*C{j}*Y'*A{j}*Y);
        end
    end
end
    fgrad = fgrad + (F+F')*Y + 2*N*(B'*B*Y-B'*b);

end

