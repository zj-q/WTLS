clear all;
clc;

%% x_ represents the parameter result
% Dx_ is its variance matrix
% i_ is iteration number 
% t_ is running time

% copy right by Zhijun Qi 7/13/2023

%% 
global y A std0 num
num=100;% number of point, A multiple of 4
%% produce the observations with errors
[std0 observation]=produce_data(num);

%% functional model
y=observation(:,2); %observation vector
A=[observation(:,1) ones(num,1)]; % coefficient matrix

%% LS method
[x_ls Dx_ls t_ls]=LS();

%% Fang2013 method
[x_wtls Dx_wtls i_wtls t_wtls]=WTLS();

%% STLS mrthod 
[x_stls Dx_stls i_stls t_stls]=STLS();

%% PEIV model in Xu2015
[x_wtls i_wtls t_wtls]=WTLS2();%% Dx_stls is not mentioned in the paper

%% IWLS in kang2022
[x_iwls Dx_iwls i_iwls t_iwls]=IWLS();

%% RTLS1 in the this paper
[x_rtls1  i_rtls1 t_rtls1]=RTLS1();

%% RTLS2 in the this paper
[x_rtls2 Dx_rtls2 i_rtls2 t_rtls2]=RTLS2();
