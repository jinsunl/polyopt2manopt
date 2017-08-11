function [ x ] = mss_delete( x,IDX )
%This function helps to delete elements in msspoly vector x
%   Input:  x -- a vector of msspoly
%           idx -- a list of index of elements that will be deleted in x
%   Output: x -- new vector of msspoly with elements deleted
IDX = sort(IDX);
for q = 1:length(IDX)
    idx = IDX(q);
    n = length(x);
    if n==1,
        x = [];
    elseif idx==1,
        x = x(2:n);
    elseif idx==n,
        x = x(1:n-1);
    else
        x = [x(1:idx-1);x(idx+1:n)];
    end
    IDX = IDX-1;
end
end
