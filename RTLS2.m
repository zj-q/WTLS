%% This function corresponds to the exact solution in Algorithm 1 in the paper

function [x_rtls2 Dx_rtls2 i_rtls2 t_rtls2]=RTLS2()
global y A std0 num

%% all the observations are divided into four groups.
num1=num/4;
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

Ea1=zeros(num1,2);
Ea2=zeros(num1,2);
Ea3=zeros(num1,2);
Ea4=zeros(num1,2);

x_2=(A'*A)\A'*y;
dx=1e5;
C1=[x_2(1)*eye(num1) -eye(num1)];

iter1=0;
%% estimate the parameter for the first data
tic
while (dx>1e-8&&iter1<100)
    N=(A1+Ea1)'*inv(C1*Q1*C1')*(A1+Ea1);
    w=(A1+Ea1)'*inv(C1*Q1*C1')*(y1+Ea1*x_2);
    x11=N\w;
    v=-Q1*C1'*inv(C1*Q1*C1')*(A1*x11-y1);
    dx=norm(x11-x_2,2);
    x_2=x11;
    Ea=[v(1:num1) zeros(num1,1)];
    C=[x_2(1)*eye(num1) -eye(num1)];
    iter1=iter1+1;
end
Qx11=inv(N);

%% estimate the parameter for the second data and previous results
y2=[y2; x11];
Q2=blkdiag(Q2,Qx11);
C2=blkdiag([x11(1)*eye(num1) -eye(num1)],eye(2));
A2=[A2;eye(2)];
Ea2=[Ea2;zeros(2,2)];
dx=10;
x_2=x11;
iter2=0;

while (dx>1e-8&&iter2<100)
    N=(A2+Ea2)'*inv(C2*Q2*C2')*(A2+Ea2);
    w=(A2+Ea2)'*inv(C2*Q2*C2')*(y2+Ea2*x_2);
    x22=N\w;
    v=-Q2*C2'*inv(C2*Q2*C2')*(A2*x22-y2);
    dx=norm(x22-x_2,2);
    x_2=x22;
    Ea2(1:num1,:)=[v(1:num1) zeros(num1,1)];     
    C2=blkdiag([x22(1)*eye(num1) -eye(num1)],eye(2));
    iter2=iter2+1;
end
Q2=diag([std0(num1+1:2*num1,1);std0(num1+1:2*num1,2)]);
Qx22=inv(N);



%%estimate the parameter for the third data and previous results
y3=[y3; x22];
Q3=blkdiag(Q3,Qx22);
C3=blkdiag([x22(1)*eye(num1) -eye(num1)],eye(2));
A3=[A3;eye(2)];
Ea3=[Ea3;zeros(2,2)];
dx=10;
x_2=x22;
iter3=0;

while (dx>1e-8&&iter3<100)
    N=(A3+Ea3)'*inv(C3*Q3*C3')*(A3+Ea3);
    w=(A3+Ea3)'*inv(C3*Q3*C3')*(y3+Ea3*x_2);
    x33=N\w;
    v=-Q3*C3'*inv(C3*Q3*C3')*(A3*x33-y3);
    dx=norm(x33-x_2,2);
    x_2=x33;
    Ea3(1:num1,:)=[v(1:num1) zeros(num1,1)];     
    C3=blkdiag([x33(1)*eye(num1) -eye(num1)],eye(2));
    iter3=iter3+1;
end
Q3=diag([std0(2*num1+1:3*num1,1);std0(2*num1+1:3*num1,2)]);
Qx33=inv(N);

%% estimate the parameter for the forth data and previous results
y4=[y4; x33];
Q4=blkdiag(Q4,Qx33);
C4=blkdiag([x33(1)*eye(num1) -eye(num1)],eye(2));
A4=[A4;eye(2)];
Ea4=[Ea4;zeros(2,2)];
dx=10;
x_2=x33;
iter4=0;

while (dx>1e-8&&iter3<100)
    N=(A4+Ea4)'*inv(C4*Q4*C4')*(A4+Ea4);
    w=(A4+Ea4)'*inv(C4*Q4*C4')*(y4+Ea4*x_2);
    x44=N\w;
    v=-Q4*C4'*inv(C4*Q4*C4')*(A4*x44-y4);
    dx=norm(x44-x_2,2);
    x_2=x44;
    Ea4(1:num1,:)=[v(1:num1) zeros(num1,1)];     
    C4=blkdiag([x44(1)*eye(num1) -eye(num1)],eye(2));
    iter4=iter4+1;
end
t_rtls2=toc;
x_rtls2=x44;
i_rtls2=iter1+iter2+iter3+iter4;
sigma=v'*inv(Q4)*v/(2*num1);
Dx_rtls2=sigma*inv(N);
end