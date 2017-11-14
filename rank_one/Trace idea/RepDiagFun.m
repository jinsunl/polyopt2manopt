function [ res ] = RepDiagFun( y,len )
%Example: if len = 3, then
% res = [ y, 0, 0 ;
%         0, y, 0 ;
%         0, 0, y ]

if size(y,2)==1,
    temp = [zeros(length(y)*len,1);y];
    temp = [y;repmat(temp,(len-1),1)];
    res = reshape(temp,[],len);
else
    temp = [y;zeros(size(y,1)*(len-1),size(y,2))];
    temp = [reshape(temp,[],1);zeros(size(y,1),1)];
    temp = repmat(temp,len,1);
    temp(end-len*size(y,1)+1:end)=[];
    res = reshape(temp,[],len*size(y,2));
end

    
end

