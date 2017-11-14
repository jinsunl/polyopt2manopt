% retraction test
clear all,close all,clc;

y = sym('y',1,'real');
Y = [1;y;y^2;y^3;y^4];
Z = [0;-y;1;-y;1]*0.05;

% Y = [1;sym('y',[4,1],'real')];
% Z = [0;sym('z',[4,1],'real')];

w = Y+Z;
for i = 1:3,
    A{i} = zeros(5,5);
    A{i}(1,i+2) = -1;
    A{i}(2,i+1) = 1;
end

%% 
yy = 2;
Y = double(subs(Y,y,yy));
Z = double(subs(Z,y,yy));
w = double(subs(w,y,yy));
alpha = sym('a',[3,1],'real');
eqn = [];
for m = 1:3,
    temp = w'*A{m}*w;
    for i = 1:3,
        temp = temp - alpha(i)*Y'*A{i}*A{m}*w - alpha(i)*w'*A{m}*A{i}'*Y;
        for j = 1:3,
            temp = temp + alpha(i)*alpha(j)*Y'*A{i}*A{m}*A{j}'*Y;
        end
    end
    eqn = [eqn;temp];
end

%%
% sol = solve(eqn,alpha,'ReturnConditions', true);
sol = vpasolve(eqn,alpha);
a1 = double(sol.a1);
a2 = double(sol.a2);
a3 = double(sol.a3);
% double(subs(eqn,alpha,[a1(1);a2(1);a3(1)]))
% double(subs(eqn,alpha,[a1(2);a2(2);a3(2)]))

alpha1 = [a1(1);a2(1);a3(1)];
alpha2 = [a1(2);a2(2);a3(2)];

res1  = w;
res2 = w;
for i = 1:3,
    res1 = res1 - A{i}'*Y*alpha1(i);
    res2 = res2 - A{i}'*Y*alpha2(i);
end
res1
res2








