function [monos,c]=p2d_decomp_noCoeff(poly)
% this function is similar to p2d_decomp except that we ignore the
% coefficient of each msspoly here. 
   [a,b,c]=decomp(poly);
   monos=msspoly(ones(size(c,2),1));
   for i=1:size(c,2)
       monos(i,1)=recomp(a,b(i,:),1);
   end
end
