%% This sub_function produces observations with errors 

function [std observation]=produce_data(num)
a1=2;% true value of slope
a2=-1;% true value of intercept

x=10*rand(num,1);
y=a1*x+a2;% true values

std1=0.5+0.5*rand(num,1);
e_error=0.1*std1.*randn(num,1);
x=e_error+x;

std2=0.5+0.5*rand(num,1);
e_error=0.1*std2.*randn(num,1);
y=e_error+y;

observation=[x y];% observation with errors
std=[std1 std2];% the standard deviation of observation
end
