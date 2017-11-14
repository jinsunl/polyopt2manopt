function [ success ] = test_auto( dim,degree,N,a,FF,randA )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    profile off
    
    x = msspoly('x',dim);
    vec = monomials(x,0:degree/2);
    temp = randn(length(vec),round(length(vec)/2));
    temp = temp*temp';
    % load temp
    f = vec'*temp*vec+8;

    % h = [1-x;x+1];
    h = [];


    %% spotless
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,lambda] = prog.newFree(1);
    prog = prog.withSOS(f-lambda);
    options = spot_sdp_default_options();
    options.verbose = 2;
    sol = prog.minimize(-lambda, @spot_mosek_sos, options);
    res_spot = sol.eval(lambda)
    pause


    %% manifold
    
    profile on
    
    if isempty(h),
            dh = 0;
        else
            dh = full(deg(h,x));
    end
    i = round(max(full(deg(f,x)),dh)/2);

    
    if randA
        F = temp;
        Mandim = nchoosek(dim+degree/2,degree/2);
        gg = floor(abs(randn(1)))+10;
        for i = 1:gg,
            A{i} = randn(Mandim);
            A{i} = 0.5*(A{i}+A{i}');
            C{i} = zeros(gg,1);
            C{i}(i) = 1;
        end
    else
%         [ F,A,C,Mandim,~] = pre_process_new( f,h,x,i,dim ); % msspoly version
%     % more pre-process needed at this moment
%         if ~isempty(A),
%             temp_idx = find(C{end}==1);
%             for i = 1:length(C),
%                 if length(C{i})>temp_idx+1,
%                     C{i}(temp_idx+2:end) = [];
%                 end
%                 A{i} = 0.5*(A{i}+A{i}');
%             end
%             b = zeros(length(A),1);
%         end
%         stackA = [];
%         A_stack = [];
%         C_stack = [];
%         for i = 1:length(A),
%             stackA = [stackA;A{i}];
%             A_stack = [A_stack;A{i}(:)'];
%             C_stack = [C_stack,C{i}];
%         end
%         CC = C_stack'*C_stack;
%         dimA = size(A{1},1);



        [ A_stack,stackA,C_stack,CC,dimA ] = pre_process_new2( dim, degree/2 ); Mandim = dimA;
    end
    

    if FF==0,
        F = 0.5*(F+F');
    else
        F = temp;F(1,1)=F(1,1)+8;
    end

    %%%%%%%%% parameters %%%%%%%%%%%%%%%
    Rank = 1;                          %
    % N =  100000000000000;                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if a==0,
        manifold = euclideanfactory(Mandim,Rank);
    else
        manifold = modefied_euclideanfactory(Mandim,Rank);
    end
    
    
    
%     %%%%%%%% VERSION 1(old) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tempA = A; tempC = C;
%     clear A; clear C
%     for i = 1:length(tempA),
%         A(:,:,i) = tempA{i};
%         C(:,:,i) = tempC{i};
%     end
%         
%     
%     % manifold = symfixedrankYYfactory(Mandim,Rank);
%     problem.M = manifold;
% 
%     % generate cost function, gradient and hessian
%     % part1: cost function
%     problem.cost = @(Y)cost( Y,F,A,C,N,b );
% 
%     % part2: gradient
%     problem.egrad = @(Y)egrad(Y,F,A,C,N,b);
% 
%     % % part3: Hessian
%     problem.ehess = @(Y,U)ehess(Y,U,F,A,C,N,b);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        
    %%%%%%%% VERSION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % manifold = symfixedrankYYfactory(Mandim,Rank);
    problem.M = manifold;

    % generate cost function, gradient and hessian
    % part1: cost function
    problem.cost = @(Y)cost_new( Y,F,A_stack,C_stack,N);

    % part2: gradient
    problem.egrad = @(Y)egrad_new(Y,F,A_stack,stackA,CC,N,dimA);

    % % part3: Hessian
    problem.ehess = @(Y,U)ehess_new(Y,U,F,A_stack,stackA,CC,N,dimA);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    problem.rhess_aug = @(Y,U)rhess_aug(Y,U,F,A_stack,stackA,N,dimA);
    problem.constrnum = size(stackA,1)/dimA;
    problem.stackA = stackA;
    problem.N = N;
    options.aug = 1;
    
    
    
    
    
    
    
    
    
    
    % % manopt
%     checkgradient(problem); 
%     pause;
    checkhessian(problem); 
%     pause;

    options.maxiter = 100000;

    [Y,cost_obj] = trustregions(problem,[],options);
    res = sign(Y(1))*Y
    cost_obj = Y'*F*Y
    cost_spot = res_spot
% 
%     if abs(cost_obj-double(res_spot))<=0.001,
%         success = 1;
%     else
%         success = 0;
%     end

    profile viewer
end



function direchess = rhess_aug(Y,U,F,A_stack,stackA,N,dimA)
    lengthA = size(stackA,1)/dimA;
    Y = Y(1:end-lengthA);
    gradC = 2*stackA'*RepDiagFun(Y,lengthA);
    M = Y*Y';
    Mtilde = A_stack*M(:);
    hessmat = [2*F+4*N*DIAGfun(Mtilde',dimA)*stackA, gradC; gradC', -1/2/N*eye(size(gradC,2))];
    direchess = hessmat * U;
    direchess(1) = 0;
end
