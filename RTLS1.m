%% This function corresponds to the exact solution in Algorithm 1 in the paper

function [x_rtls1  i_rtls1 t_rtls1]=RTLS1()
global y A std0 num
num1=num/4;% all the observations are divided into four groups.
Ea1=zeros(num1,2);
Ea2=zeros(num1,2);
Ea3=zeros(num1,2);
Ea4=zeros(num1,2);

y1=y(1:num1);
y2=y(num1+1:2*num1);
y3=y(num1*2+1:3*num1);
y4=y(3*num1+1:4*num1);

A1=A(1:num1,:);
A2=A(num1+1:2*num1,:);
A3=A(num1*2+1:3*num1,:);
A4=A(3*num1+1:4*num1,:);

Q1=diag([std0(1:num1,1);std0(1:num1,2)]);
Q2=diag([std0(num1+1:2*num1,1);std0(num1+1:2*num1,2)]);
Q3=diag([std0(2*num1+1:3*num1,1);std0(2*num1+1:3*num1,2)]);
Q4=diag([std0(3*num1+1:4*num1,1);std0(3*num1+1:4*num1,2)]);

x_1=(A'*A)\A'*y;
C1=[x_1(1)*eye(num1) -eye(num1)];
dx=100;
iter1=0;

%% Begin the iteration
tic
while (dx>1e-8&&iter1<100)
    N=(A1+Ea1)'*inv(C1*Q1*C1')*(A1+Ea1);
    w=(A1+Ea1)'*inv(C1*Q1*C1')*(y1+Ea1*x_1);
    x1=N\w;
    Qx1=inv(N);

    % begin first recursion
    J2=inv((A2+Ea2)'*inv(C1*Q2*C1')*A2+inv(Qx1))*(A2+Ea2)'*inv(C1*Q2*C1');
    v2=y2-A2*x1;
    x2=x1+J2*v2;
    Qx2=inv(inv(Qx1)+(A2+Ea2)'*inv(C1*Q2*C1')*(A2+Ea2));

    % begin second recursion
    J3=inv((A3+Ea3)'*inv(C1*Q3*C1')*A3+inv(Qx2))*(A3+Ea3)'*inv(C1*Q3*C1');
    v3=y3-A3*x2;
    x3=x2+J3*v3;
    Qx3=inv(inv(Qx2)+(A3+Ea3)'*inv(C1*Q3*C1')*(A3+Ea3));

    % begin forth recursion
    J4=inv((A4+Ea4)'*inv(C1*Q4*C1')*A4+inv(Qx3))*(A4+Ea4)'*inv(C1*Q4*C1');
    v4=y4-A4*x3;
    x4=x3+J4*v4;
    Qx4=inv(inv(Qx3)+(A4+Ea4)'*inv(C1*Q4*C1')*(A4+Ea4));

    %% renewed the vector and matrix
    v11=-Q1*C1'*inv(C1*Q1*C1')*(A1*x4-y1);
    v22=-Q2*C1'*inv(C1*Q2*C1')*(A2*x4-y2);
    v33=-Q3*C1'*inv(C1*Q3*C1')*(A3*x4-y3);
    v44=-Q4*C1'*inv(C1*Q4*C1')*(A4*x4-y4);
    
    Ea1=[v11(1:num1) zeros(num1,1)];
    Ea2=[v22(1:num1) zeros(num1,1)];
    Ea3=[v33(1:num1) zeros(num1,1)];
    Ea4=[v44(1:num1) zeros(num1,1)];

    dx=norm(x4-x_1,2);
    x_1=x4;
    C1=[x4(1)*eye(num1) -eye(num1)];

    iter1=iter1+1;
end
i_rtls1=iter1;
x_rtls1=x4;
t_rtls1=toc;
end
