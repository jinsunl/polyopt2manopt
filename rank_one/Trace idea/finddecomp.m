function [ i,j ] = finddecomp( M,q )
%this function helps to find the decomposition of q in M s.t.
%                      q = M(i)*M(j)
%where M is a vector of msspoly, and q is a single msspoly.                
%if no decomposition is found, then i = j = 0

%     [x,p,m] = decomp(q);
    n = length(M);
    i = 0;
    j = 0;
    found = 0;
    for ii = n:-1:2,
        for jj = ii:-1:2,
            if isequal(q,M(ii)*M(jj)),
                i = ii;
                j = jj;
                found = 1;
                break;
            end
        end
        if found,
            break;
        end
    end

end

