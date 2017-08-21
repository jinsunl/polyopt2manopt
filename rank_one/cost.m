function [ f ] = cost( Y,F,A,C,N,B,b )
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates
    
    temp = 0;
if ~isempty(A),
    for i = 1:length(A),
       temp = temp + C{i}*Y'*A{i}*Y; 
    end
end
    f = Y'*F*Y+N*norm(temp+B*Y-b,'fro')^2;


end

