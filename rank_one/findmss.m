function [ ii,jj ] = findmss( M,x )
%this function helps to find the index of x in a matrix M. If x cannot be
%found in M, then ii = jj = 0;

[row,col] = size(M);
ii = 0;
jj = 0;
found = 0;
if row==col && isequal(M,M'), % M is symmetric
    for i = 1:col,
        for j = i:row
            if isequal(x,M(j,i)),
               ii = i;
               jj = j;
               found = 1;
               break;
            end
        end
        if found,
            break;
        end
    end

else    % M is not symmetric
    for i = 1:col
        for j = 1:row
            if isequal(x,M(j,i)),
               ii = i;
               jj = j;
               found = 1;
               break;
            end
        end
        if found,
            break;
        end
    end
end

end

