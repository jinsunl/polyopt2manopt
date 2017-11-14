function [ f ] = cost( Y,F,A,C,N,b )
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates
    
    temp = 0;
    if ~isempty(A),
        for i = 1:size(A,3),
           temp = temp + C(:,:,i)*trace(A(:,:,i)*Y*Y'); 
        end
    end
    f = trace(F*Y*Y')+N*norm(temp-b,'fro')^2;


end

