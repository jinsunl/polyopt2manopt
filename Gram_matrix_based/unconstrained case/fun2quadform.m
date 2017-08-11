function [ MM,Basis,Mono ] = fun2quadform( F, var)
%This function transfer polynomials into matrix representation
%   Input:  F -- a vector of polynomials
%           var -- a vector of indeterminate variables in F
%   Output: M -- cells of matrices such that F(i)=basis{i}'*M{i}*basis{i}
%           Basis -- cells of vectors of monomials of var, which is needed
%                    to represent F
%           Mono -- cells of vectors of monomials in F(i)


%     [x,p,M] = decomp(f);
%     cont = nchoosek(length(x)+round(d/2),length(x));

for q = 1:length(F),
    f = F(q);
    nvar = length(var);
    d = full(deg(f,var));
    basis = monomials(var,0:round(d/2));
    Q = msspoly(zeros(length(basis)));
    mono = p2d_decomp(f);
    monos.real = mono;
    monos.base = mono;
    monos.para = mono;
    for i = 1:length(mono),
       temp = var;
       [x,p,M] = decomp(mono(i));
       list = 1:length(x);
       listvar = [];
       for ii = 1:length(x),
           for jj = 1:length(temp),
              if isequal(x(ii),temp(jj)),
                  listvar = [listvar,ii];
                  temp = mss_delete(temp,jj);
                  break;
              end
           end
       end
       listfree = setdiff(list,listvar);
       
       monos.base(i) = msspoly(1);
       for ii = 1:length(listvar),
           monos.base(i) = monos.base(i)*x(listvar(ii))^p(1,listvar(ii));
       end
       
       monos.para(i) = msspoly(M);
       for ii = 1:length(listfree),
           monos.para(i) = monos.para(i)*x(listfree(ii))^p(1,listfree(ii));
       end
    end
%     monos

    list = 1:length(mono);
    for i = 1:length(basis),
       k = length(list);
       if k==0,
           break;
       end
       while k>0,
          if isequal(basis(i)*basis(i), monos.base(list(k)))
              Q(i,i) = Q(i,i) + monos.para(list(k));
              list(k) = [];
          end
          k = k-1;
       end
       if length(list)==0,
           break;
       end
    end
    
    for i = 1:length(basis)-1,
        for j = i+1:length(basis),
            k = length(list);
            if k==0,
                break;
            end
            while k>0
                if isequal(basis(i)*basis(j), monos.base(list(k)))
                        Q(i,j) = Q(i,j)+monos.para(list(k))/2;
                        Q(j,i) = Q(i,j);
                    list(k) = [];
                end
                k = k-1;    
            end
        end
        if length(list)==0,
            break;
        end
    end
    MM{q} = Q;
    Basis{q} = basis;
    Mono{q} = mono;
    
    
end
end

