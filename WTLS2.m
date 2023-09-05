%% Alternative formulae for parameter estimation in partial errorsin-variables model,Journal of geodesy

function [x_wtls2 i_wtls2 t_wtls2]=WTLS2()
global y A std0 num
Qa=diag(std0(:,1));% cofactor matrix for matrix A
Qy=diag(std0(:,2));% cofactor matrix for vector y

%% fixing vector and matrix
h=[zeros(num,1);A(:,2)];
B=[eye(num);zeros(num,num)];

a = A(:,1);
aa = A(:,1);
x_wtls2 = zeros(2,1);
iter = 0;

tic
for k = 1:1000
    iter=iter+1;
    AA=reshape(h+B*a,num,2);
    xx=(AA'*(Qy\AA))\(AA'*(Qy\y));
    S=kron(xx',eye(num))*B;
    E = Qy + S*Qa*S';
    a = aa + Qa*S'*(E\(y-A*xx));
    if norm(xx-x_wtls2)<10^(-8)
        break;
    end
     x_wtls2=xx;
end
t_wtls2=toc;
i_wtls2=iter;
%% accuracy assessment
ea=aa-a;
ey=y-AA*xx;
sigma=(ea'*inv(Qa)*ea+ey'*inv(Qy)*ey)/(num*2-2);
end
