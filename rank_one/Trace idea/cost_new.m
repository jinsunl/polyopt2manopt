function [ f ] = cost_new( Y,F,A_stack,C_stack,N)
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates
    
    temp = 0;
    M = Y*Y';
    if ~isempty(A_stack),
        temp = C_stack*(A_stack*M(:));
    end
    f = Y'*F*Y + N*norm(temp,'fro')^2;


end

