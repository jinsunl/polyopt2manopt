function [ base,para ] = para_mono( mono,var )
%This function helps to transfer vector of monomials w.r.t. all variales
%into vector of monomials w.r.t. var only.
%   Example:
%   mono = [           (6)  ]     var = [   x1  ]
%          [       (-2)*v1  ]
%          [        (4)*x1  ]
%          [      (6)*x1^2  ]
%          [      (4)*x1^3  ]
%          [          x1^4  ]
%          [    (-2)*v2*x1  ]
%          [  (-2)*v3*x1^2  ]
%          [      -v4*x1^2  ]
%          [  (-2)*v5*x1^3  ]
%          [      -v6*x1^4  ]
%   then
%   base = [   (1)  ]            para = [     (6)+(-2)*v1  ]
%          [    x1  ]                   [     (4)+(-2)*v2  ]
%          [  x1^2  ]                   [  (6)+(-2)*v3-v4  ]
%          [  x1^3  ]                   [     (4)+(-2)*v5  ]
%          [  x1^4  ]                   [          (1)-v6  ]
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
    
    base = monos.base(1);
    para = monos.para(1);
    for i = 1:length(monos.base),
        if i==1,
            continue;
        end
        flag = 1;
        for j = 1:length(base),
           if isequal(base(j),monos.base(i)),
               para(j) = para(j)+monos.para(i);
               flag = 0;
               break;
           end
        end
        if flag==1,
           base(end+1,1) = monos.base(i);
           para(end+1,1) = monos.para(i);
        end
        
        
    end
    
end

