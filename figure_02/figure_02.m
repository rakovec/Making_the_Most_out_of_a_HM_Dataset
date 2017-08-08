%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the conditional and unconditional model output densities for FUSE-016 model
%
% written by Xuefei LU
% 2016 Dec
% email: xuefei.lu@unibocconi.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% import data
currentpath = cd('..');
addpath('./data_input')
addpath('./figure_02')

modelnr = 016;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);

%% Visualization of estimated conditional and unconditional model output densities
figure
vnames = labels(IND);
[d,s,p]=deltafast2(x(:,IND),y,50,'cond',vnames);