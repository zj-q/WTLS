function [x_ls Dx_ls t_ls]=LS()
global y A std0 num
Q1=diag(std0(:,2));
P=inv(Q1);
tic
x_ls=(A'*P*A)\(A'*P*y);
t_ls=toc;
v=y-A*x_ls;
sigma0=v'*P*v/(num-2);
Dx_ls=sigma0*inv(A'*P*A);
end