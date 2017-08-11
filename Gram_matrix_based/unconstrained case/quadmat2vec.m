function [ M_tilde,M_const ] = quadmat2vec( M,var )
%This function vectorize matrix M{i} containing var into the form as:
%          vec(M{i}) = M_tilde{i}*var+M_const{i}
%where M_tilde is a matrix and M_const is a vector.
%   Input:  M -- cells of matrices 
%           var -- vector of variables that lie inside M
%   OUtput: M_tilde -- cells of matrices
%           M_const -- cells of vectors
for q = 1:length(M),
    mat = M{q};
    fun = fn(mat,var);
    x = sym('x',size(var),'real');
    mat = fun(x);
    
    n = size(mat,1);
    x_aux = sym('v',[n*(1+n)/2,1],'real');
    idx = 1;
    mat0 = sym(zeros(n));
    for i = 1:n,
        for j = i:n,
            mat0(i,j) = x_aux(idx);
            idx = idx+1;
        end
    end
    mat0 = mat0 + triu(mat0,1)';
    mat = mat+mat0;
    x = [x;x_aux];
    
    mat = reshape(mat,size(mat,1)*size(mat,2),1);
    M_tilde{q} = double(jacobian(mat,x));
    M_const{q} = double(mat-M_tilde{q}*x);
    
    
    
end

end

