clear all,close all,clc

%% projection test
dim = 3;
degr = 2;
y = msspoly('y',dim);
Yvec = monomials(y,0:degr);
Mdim = length(Yvec);

temp = zeros(Mdim); temp(2,2)=1;temp(1,3)=-1;
A{1} = temp;
temp = zeros(Mdim); temp(2,4)=1;temp(1,5)=-1;
A{2} = temp;
temp = zeros(Mdim); temp(4,4)=1;temp(1,6)=-1;
A{3} = temp;
temp = zeros(Mdim); temp(2,7)=1;temp(1,8)=-1;
A{4} = temp;
temp = zeros(Mdim); temp(4,7)=1;temp(1,9)=-1;
A{5} = temp;
temp = zeros(Mdim); temp(7,7)=1;temp(1,10)=-1;
A{6} = temp;
% temp = zeros(Mdim); temp(2,11)=1;temp(1,12)=-1;
% A{7} = temp;
% temp = zeros(Mdim); temp(4,11)=1;temp(1,13)=-1;
% A{8} = temp;
% temp = zeros(Mdim); temp(7,11)=1;temp(1,14)=-1;
% A{9} = temp;
% temp = zeros(Mdim); temp(11,11)=1;temp(1,15)=-1;
% A{10} = temp;
% temp = zeros(Mdim); temp(2,16)=1;temp(1,17)=-1;
% A{11} = temp;
% temp = zeros(Mdim); temp(4,16)=1;temp(1,18)=-1;
% A{12} = temp;
% temp = zeros(Mdim); temp(7,16)=1;temp(1,19)=-1;
% A{13} = temp;
% temp = zeros(Mdim); temp(11,16)=1;temp(1,20)=-1;
% A{14} = temp;
% temp = zeros(Mdim); temp(16,16)=1;temp(1,21)=-1;
% A{15} = temp;


for i = 1:length(A),
    A{i} = sym_mat(A{i});
end




Y = double(subs(Yvec,y, randn(dim,1)));
z = randn(Mdim,1);
e1 = zeros(Mdim,1);e1(1) = 1;
temp = [];
for j = 1:length(A),
    for i = 1:length(A),
        AA(j,i) = Y'*A{j}*A{i}*Y;
    end
    AA(j,i+1) = Y'*A{j}*e1;
    B(j,1) = Y'*A{j}*z;
    temp = [temp,e1'*A{j}*Y];
end
AA = [AA;temp,1];
B = [B;e1'*z];
rank(AA)
alpha = AA\B;
res = z-alpha(end)*e1;
for i = 1:length(A),
    res = res-alpha(i)*A{i}*Y;
end

%%
% Z = res;
% Ybar = Y+Z;
% j = 9;
% temp = zeros(17);
% temp(1,1) = Ybar'*A{j}*Ybar;
% for i = 1:15,
%     temp(1,i+1) = Ybar'*A{i}*A{j}*Ybar;
%     temp(i+1,1) = temp(1,i+1);
%     temp(17,i+1) = Ybar'*A{i}*A{j}*e1;
%     temp(i+1,17) = temp(17,i+1);
% end
% temp(1,17) = e1'*A{j}*Ybar;
% temp(17,1) = temp(1,17);
% temp(17,17) = e1'*A{j}*e1;
% for m = 1:15,
%     for n = 1:15,
%         temp(m+1,n+1) = Ybar'*A{m}*A{j}*A{n}*Ybar;
%     end
% end
% eig(temp)
% spy(temp)





%% test retraction (works, but not fast)
Z = res;
Ybar = Y+Z;
alpha = sym('a',[length(A),1],'real');
b = sym('b','real');
func = sym(zeros(length(A)+1,1));
numA= length(A);

for i = 1:numA,
    func(1) = func(1) - alpha(i)*e1'*A{i}*Ybar;
end
func(1) = func(1)+e1'*Ybar-b==1;

for j = 1:numA,
    func(j+1) = Ybar'*A{j}*Ybar-2*b*e1'*A{j}*Ybar+b^2*e1'*A{j}*e1;
    for m = 1:numA,
        for n = 1:numA,
            func(j+1) = func(j+1)+alpha(m)*alpha(n)*Ybar'*A{m}*A{j}*A{n}*Ybar;
        end
        func(j+1) = func(j+1)-2*alpha(m)*Ybar'*A{m}*A{j}*Ybar+2*alpha(m)*b*Ybar'*A{m}*A{j}*e1;
    end
    func(j+1) = func(j+1)==0;
end

solution = vpasolve(func);  % not unique solution
% solution = solve(func);   % not working
alpha = double([solution.a1,solution.a2,solution.a3,solution.a4,solution.a5,solution.a6]);
b = double(solution.b);

Ynew_list = ones(1,length(b));
for j = 1:length(b),
    Ynew{j} = Ybar - b(j)*e1;
    for i = 1:numA,
        Ynew{j} = Ynew{j} - alpha(j,i)*A{i}*Ybar;
    end
    [Y Ybar Ynew{j} Yvec]
    if ~IsOnManifold(A,Ynew{j}),
        Ynew_list(j)=0;
    end
end
Ynew_list
temp = [alpha b];
temp = sum(abs(temp).^2,2).^(1/2);
Ynew_real = Ynew{find(temp==min(temp))}



%% help function
function mat = sym_mat(mat)
    mat = 0.5*(mat+mat');
end

    
    
function res = IsOnManifold(A,Y)
   res = 1;
   for i = 1:length(A),
      if abs(Y'*A{i}*Y)>1e-8,
          res = 0;
          break
      end
   end
end






