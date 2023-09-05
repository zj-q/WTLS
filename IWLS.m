%% this function is iterative weight least squares

function [x_iwls Dx_iwls i_iwls t_iwls]=IWLS()
global y A std0 num
Qa=diag([std0(:,1);zeros(num,1)]);% cofactor matrix for matrix A
Qy=diag(std0(:,2));% cofactor matrix for vector y

x_first=(A'*A)\A'*y;
dx=10;
iter=0;
tic
while (dx>1e-8&&iter<1000)
    XX=kron(x_first,eye(num));
    M=inv(Qy+XX'*Qa*XX);
    x_iwls=(A'*M*A)\A'*M*y;
    dx=norm(x_iwls-x_first);
    x_first=x_iwls;
    iter=iter+1;
end
t_iwls=toc;
i_iwls=iter;
v=y-A*x_iwls;
sigma=v'*M*v/(num-2);
Dx_iwls=sigma*inv(A'*M*A);
end

