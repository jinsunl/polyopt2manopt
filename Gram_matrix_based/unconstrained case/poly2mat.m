function [ M,Const,b_aux, Xdim ] = poly2mat( F, var, vlambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for q = 1:length(F),
    f = F(q);
    d = full(deg(f,var));
    basis = monomials(var,0:round(d/2));
    Xdim{q} = length(basis);
    n = length(basis);
    vG = msspoly('v',n*(n+1)/2);
    v = [vlambda;vG]; % 'decision variable' when spotless is in use
    G = msspoly(zeros(n));% Gram matrix which represent f
    idx = 1;
    for i=1:n,
        for j = 1:n,
            if j>=i,
                G(i,j) = vG(idx);
                idx = idx+1;
            else
                G(i,j) = G(j,i);
            end
        end
    end
    
    % construct M1 s.t. M1*v=vec(YY') where Y is our decision variable on the manifold. 
    f_gram = basis'*G*basis;
    m = sym('m',[length(v),1],'real');
    fun = fn(G,v);
    Gsym = fun(m);
    Gsym = reshape(Gsym,size(Gsym,1)^2,1);
    temp = [];
    for i = 1:n,
        temp = [temp,(i-1)*n+[1:n]'];
    end
    list = tril(temp);
    list = reshape(list,n^2,1);
    del_list = find(list==0);
    list(del_list)=[];
    M1 = double(jacobian(Gsym(list),m));
        
    
    
     
    
    
    
    
    % construct M2,Const s.t. M2*v+Const = 0. This represent f-f_gram=0,
    % which means G is the gram matrix of f. 
    f_zero = f-f_gram;
%     mono_fzero = p2d_decomp(f_zero);
%     [base,para] = para_mono(mono_fzero,var);
    [base,para] = mss2para_mono(f_zero,v);
    
    fun = fn(para,v);
    para = fun(m);
    M2 = double(jacobian(para,m));
    Const{q} = [zeros(n*(n+1)/2,1);double(para-M2*m)];
    
    % construct M and b_aux s.t. b_aus*vec(YY') = [vec(YY');0]
    b_aux{q} = [eye(n^2);zeros(length(para),n^2)];
    b_aux{q}(del_list,:) = [];
    M{q} = [M1;M2];
    
    
    
    
end
end

