function [ F,A,C,Mandim,B,b ] = pre_process( f,h,x,i )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    


    df = full(deg(f,x));
    if isempty(h),
        dh = 0;
    else
        dh = full(deg(h,x));
    end
    d = max(df,dh);
    n = d/2;
    nh = length(h);
    if mod(df,2) == 1;
        error('degree of f should be even.');
    end
    if i<n,
       error('i cannot be smaller than half of degree of f and h.') 
    end
    
    mono = monomials(x,0:i);
    zero = zeros(1,length(mono)+nh);
    M = mono*mono'; % moment matrix in msspoly
    F = 0;
    [monof,ch] = p2d_decomp_noCoeff(f);
    for ii = 1:length(monof),
        [idx1,idx2] = findmss(M,monof(ii));
        temp1 = zero;temp1(idx1) = 1;
        temp2 = zero;temp2(idx2) = 1;% thus monof(ii) = mono'*temp1'*temp2*mono 
        F = F+temp1'*temp2*ch(ii);
    end

%   auxilary parameters to allocate the following constraints as a vector. 
    count = 1;
    c = zeros(length(mono)-length(x)+nh+1,1);


% represent h using elements in M:
%   suppose Y is an element on rankfixYYfactory(manifold), then we hope to
%   generate the following constraint: Y'A{i}Y=0
%   we can do this because any element in M&h can be represented as Y'AY,
%   where Y can be think as [mono;sqrt(h)].
if nh>0,    
    for ii = 1:nh,
        temp3 = zero; temp3(end-nh+ii) = 1;
        A{ii} =  - temp3'*temp3;
        [monoh,ch] = p2d_decomp_noCoeff(h(ii));
        for iii = 1:length(monoh),
            [idx1,idx2] = findmss(M,monoh(iii));
            temp1 = zero;temp1(idx1) = 1;
            temp2 = zero;temp2(idx2) = 1;  % thus monoh(iii) = mono'*temp1'*temp2*mono 
            A{ii} = A{ii} + temp1'*temp2*ch(iii);
        end
        C{count} = c;
        C{count}(count) = 1;
        count = count+1;
    end
end    
% construct relations among elements of mono in Y. Again, we represent
% these relationships as Y'A{i}Y = 0;
    tempmono = [];
    for ii = 1:length(mono),
        if deg(mono(ii),x)<2,
            tempmono = [tempmono;mono(ii)];
            continue;
        end
        [idx1,idx2] = finddecomp(tempmono,mono(ii));
        tempmono = [tempmono;mono(ii)];
        temp1 = zero;temp1(idx1) = 1;
        temp2 = zero;temp2(idx2) = 1; % thus mono(ii) = mono'*temp1'*temp2*mono 
        temp3 = zero;temp3(1) = 1;
        temp4 = zero;temp4(ii) = 1; % thus mono(ii) = mono'*temp3'*temp4*mono 
        A{count} = temp1'*temp2 - temp3'*temp4;
        C{count} = c;
        C{count}(count) = 1;
        count = count+1;
    end
    
    Mandim = length(mono)+nh;
    
    if ~exist('A','var'),
        A = [];
    end
    if ~exist('C','var'),
        C = [];
    end
    
    % for constraint: Y(1) = 1
    b = c; b(end) = 1; zero(1) = 1;
    B = b * zero;
    
end

