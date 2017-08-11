function [ base,para ] = mss2para_mono( f,decvar )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [a,b,c] = decomp(f,decvar);
    para = c';
    for i=1:size(c,2)
       base(i,1)=recomp(a,b(i,:),1);
    end
    
    
    

end

