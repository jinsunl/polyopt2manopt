function [ A_stack,stackA,C_stack,CC,dimA ] = pre_process_new2( dim, degree )
%this function helps to generate magtic A matrices and corresbonding C
%matrices. A's and C's are then stacked into the outputs. 
%   dim: # of indeterminates in monomial vector (or say moment vector, like [1 x x^2 x^3 ...]')
%   degree: max degree of elements in moment vector.
%   NOTICE: no semi-algebriac constraints are considered at this time
%           degree is no greater than 35. (this 35 comes from the limitation of matlabfunction 'dec2base', and we need degree+1<=36)
    
    dimA = nchoosek(dim+degree,dim); % ideally A is a square matrix, so its size should be dimA*dimA, and dimA = the lenth of the moment vector.
    lengthA = dimA-1-dim; % this represents how many A matrices we'll have. '-1' means means no A is generated based on the '1' element in moment vector, and '-dim' means no A is generated based on the first degree terms in the moment vector. 
    
    if lengthA==0, % this means we dont have any monomials with degree>1, thus no need to generate magic A matrices. 
        A_stack = [];
        stackA = [];
        C_stack = [];
        CC = [];
        return
    end
    
% generate C (or say C_stack and CC). C_stack and CC will be used in cost, egrad and ehess.
    C_stack = eye(lengthA); % C_stack = [C1,C2,...,Cn]=eye(n) where n = lengthA. Basically we'll have lengthA polynomial-type cconstraints for the moment vector, and Ci is to stack those constraints into a single column vector. Thus C_stack happens to be identity.
    CC = eye(lengthA); % CC = C_stack'*C_stack, and happens to be identity as well.
    
% generate A (or say A_stack and stackA).
    % generate fake moment vector
    moment_vec = zeros(dimA,1); % at the end of the day, this vector will only contained decimal numbers that represent valid monomials. say we have monomials x^3y^2z, then it can be represended as '321' in the base of (degree+1). Convert '321' into decimal, and the result will be saved in moment_vec
    ptr1 = degree*(degree+1)^(dim-1);
    moment_vec_based = zeros(dimA,dim);  % this vector contains the (degree+1)-based number of valid monomials. In this way, ith row of this vector represent the degrees of the ith monomial w.r.t. each indeterminate. Again, say we have monomials x^3y^2z, then it will be saved as [3 2 1] as a row of moment_vec_based.
    ptr2 = dimA;
    while (ptr1>=0),
        temp = ptr1;
        temp_based = dec2base(temp,degree+1,dim); % change dec to (degree+1)-based value
        Alist = find(temp_based>64); temp_based(Alist) = temp_based(Alist)-65+10;
        numlist = find(temp_based>47);temp_based(numlist) = temp_based(numlist)-48;
        if sum(double(temp_based))<=degree,
            moment_vec_based(ptr2,:) = double(temp_based);
            moment_vec(ptr2) = ptr1;
            ptr2 = ptr2-1;
        end
        ptr1 = ptr1-1;
    end
    
    % generate A_stack and stackA. A{.} represents our cell version of
    % magic A matrices.
    A_stack = zeros(lengthA,dimA^2); % A_stack = [A{1}(:)';A{2}(:)';...;A{lengthA}(:)']
    stackA = zeros(lengthA*dimA,dimA); % stackA = [A{1};A{2};...;A{lengthA}]
    j = 1; % counter for magic matrix A's. I need this for stackA and stack_A
    for i = 1:dimA,
       % Let's try to decomp each monomials in moment_vec
       monomial = moment_vec_based(i,:);
       if sum(monomial)<=1, % if this happens, this means the monomial is with degree no greater than 1, so that no A{.} is needed.
           continue
       end
       
       if sum(find(monomial>1))==0,
            decomp1 = 0*monomial;
            decomp1(find(monomial==1,1)) = 1;
            decomp2 = monomial-decomp1;
       else
            decomp1 = floor(monomial/2); % we want to have monomial = decomp1*decomp2. Again, monomials decomp1 and decomp2 is represented by its corresponding degree vector. 
            decomp2 = monomial-decomp1;
       end
       
       % Now, find the index for decomp1 and decomp2 in moment_vec
       % first, we have to transfer decomp1 and decomp2 into (degree+1)based representation
       decomp1 = char(decomp1); Alist = find(decomp1>9);numlist = find(decomp1<10);
       decomp1(Alist) = decomp1(Alist) + 55; decomp1(numlist) = decomp1(numlist) + 48;
       decomp1 = base2dec(decomp1,(degree+1)); IDX1 = double(moment_vec==decomp1);
       
       decomp2 = char(decomp2); Alist = find(decomp2>9);numlist = find(decomp2<10);
       decomp2(Alist) = decomp2(Alist) + 55; decomp2(numlist) = decomp2(numlist) + 48;
       decomp2 = base2dec(decomp2,(degree+1)); IDX2 = double(moment_vec==decomp2);

       % finally, let's generate A and other stuff
       Aj = IDX1*IDX2'; Aj(1,i) = -1; Aj = 0.5*(Aj+Aj');  % Aj is our ith magic A matrix 
       A_stack(j,:) = Aj(:)';
       stackA((j-1)*dimA+1:j*dimA,:) = Aj;
       j = j+1;
    end
    
    
junk = 1;
end

