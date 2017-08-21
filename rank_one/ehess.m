function [ fhess ] = ehess(Y,U,F,A,C,N, B, b)
%This function generates the Euclidian Hessian of the cost function. (user
%design)
%   Input: user indecates
    
    fhess = 0;
if ~isempty(A),
    n = length(A);
    for i = 1:n,
        fhess = fhess + A{i}*U*Y'*B'*C{i} + B'*C{i}*U'*A{i}*Y + A{i}'*U*C{i}'*B*Y + A{i}'*Y*C{i}'*B*U+...
                        A{i}*Y*U'*B'*C{i} + B'*C{i}*Y'*A{i}*U - A{i}*U*b'*C{i} - A{i}'*U*C{i}'*b;
        for j = 1:n,
            fhess = fhess + A{i}*U*Y'*A{j}*Y*C{j}'*C{i}+...
                            A{i}*Y*Y'*A{j}*U*C{j}'*C{i}+...
                            A{j}*Y*C{j}'*C{i}*U'*A{i}*Y;
            fhess = fhess + A{i}'*U*C{i}'*C{j}*Y'*A{j}*Y+...
                            A{i}'*Y*C{i}'*C{j}*Y'*A{j}*U+...
                            A{j}*Y*U'*A{i}'*Y*C{i}'*C{j};
        end
    end
end
    fhess = 2*N*fhess + (F+F')*U+2*N*B'*B*U;
end

