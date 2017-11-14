function [ A ] = DIAGfun( a, dim )
%Example: input a=[1 3 4 5 2], dim=x
%         output A = [1eye(x),3eye(x),4eye(x),5eye(x),2eye(x)]
    temp = [zeros(dim,length(a));a];
    temp = repmat(temp,[dim-1,1]);
    a = [a;temp];
    A = reshape(a,dim,[]);

end

