%% STLS function, similar with WTLS
%% "Comparison of Structured and Weighted Total LeastSquares Adjustment Methods for Linearly Structured Errors-in-Variables Models
function [x_stls Dx_stls i_stls t_stls]=STLS()
global y A std0 num

%% statistical model
Q=diag([std0(:,1);std0(:,2)]);

iter=0;
Ea=zeros(num,2);
x_stls=(A'*A)\A'*y;%% inticial guess by LS
dx=1e5;
C=[x_stls(1)*eye(num) -eye(num)];%% characteristic matrix 
tic
%% Iteration begins
while (dx>1e-8&&iter<100)
    N=(A+Ea)'*inv(C*Q*C')*A;
    w=(A+Ea)'*inv(C*Q*C')*y;
    x=N\w;
    v=-Q*C'*inv(C*Q*C')*(A*x-y);
    dx=norm(x-x_stls,2);
    x_stls=x;
    Ea=[v(1:num) zeros(num,1)];
    C=[x_stls(1)*eye(num) -eye(num)];
    iter=iter+1;
end
i_stls=iter;
t_stls=toc;

sigma=v'*inv(Q)*v/(num*2-2);
Dx_stls=sigma*inv((A+Ea)'*inv(C*Q*C')*(A+Ea));
end