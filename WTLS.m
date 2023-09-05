%% Weighted total least squares
%% Weighted total least squares: necessary and sufficient conditions, fixed and random parameter, Journal of Geodesy
function [x_wtls Dx_wtls i_wtls t_wtls]=WTLS()
global y A std0 num

%% statistical model
Q1=diag([std0(:,1);std0(:,2)]);
K=zeros(3*num,2*num);
K(1:num,1:num)=eye(num);
K(2*num+1:3*num,num+1:2*num)=eye(num);
Q=K*Q1*K';

iter=0;% iteration index
Ea=zeros(num,2); % error matrix
x_wtls=(A'*A)\A'*y;% inticial guess
dx=1e5;
B=[kron(x_wtls',eye(num)) -eye(num)];

%% begin iteration
tic
while (dx>1e-8&&iter<100)
    N=(A+Ea)'*inv(B*Q*B')*A;
    w=(A+Ea)'*inv(B*Q*B')*y;
    x=N\w;
    v=-Q*B'*inv(B*Q*B')*(A*x-y);
    dx=norm(x-x_wtls,2);
    x_wtls=x;
    Ea=[v(1:num) v(num+1:2*num)];
    B=[kron(x_wtls',eye(num)) -eye(num)];
    iter=iter+1;
end
i_wtls=iter;
t_wtls=toc;% time for iteration

%% accuracy assessment
sigma=v'*pinv(Q)*v/(num*2-2);
Dx_wtls=sigma*inv((A+Ea)'*inv(B*Q*B')*(A+Ea));
end
