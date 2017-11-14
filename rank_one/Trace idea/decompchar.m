function [ decomp1,decomp2 ] = decompchar( str )
%this function help to achieve str=decomp1+decomp2 in the sense of numeric
%sum.
%   Detailed explanation goes here
   
    if sum(find(str>'1'))==0,
        decomp1 = char(0*str+'0');
        decomp1(find(str=='1',1))='1';
        decomp2 = char(str - decomp1+'0');
    else
        Alist = find(str>64); str(Alist) = str(Alist)-55;
        numlist = find(str>47); str(numlist) = str(numlist)-'0';
        decomp1 = char(floor(str/2)); decomp2 = char(str-decomp1);
        
        numlist = find(decomp1<10); Alist = find(decomp1>9);
        decomp1(numlist) = decomp1(numlist) + '0';
        decomp1(Alist) = decomp1(Alist) + 55;
        
        numlist = find(decomp2<10); Alist = find(decomp2>9);
        decomp2(numlist) = decomp2(numlist) + '0';
        decomp2(Alist) = decomp2(Alist) + 55;
 

        
    end
    
end

