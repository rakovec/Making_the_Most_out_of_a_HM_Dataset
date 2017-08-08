%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the Boxplots for the bootstrap estimates for FUSE-170 model
%Four sensitivity measures are considered:
% First order Sobol, Borgonovo delta, Kuiper index, median DELSA
%
% written by Xuefei LU
% 2016 Dec
% email: xuefei.lu@unibocconi.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% import data
currentpath = cd('..');
addpath('./data_input')
addpath('./figure_S1')

modelnr = 170;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);
delsarmse = csvread(['fuse',num2str(modelnr,'%03i'),'_delsa_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);

%% Boxplots for the bootstrap estimates
rng(1234)
nboot = 500;
E = zeros(nboot,k);
delta = zeros(nboot,k); 
T = zeros(nboot,k);
Delsaboot = zeros(nboot,k);

for i = 1:k
     [theta_hat_hat,theta_hat]=strapest(nboot,@bootbetaKS3eta,x(:,i),y,50);
     E(:,i) = theta_hat_hat;
end

for i = 1:k
     [theta_hat_hat,theta_hat]=strapest(nboot,@bootbetaKS3d,x(:,i),y,50);
     delta(:,i) = theta_hat_hat;
end

for i = 1:k
     [theta_hat_hat,theta_hat]=strapest(nboot,@bootbetaKS3k,x(:,i),y,50);
     T(:,i) = theta_hat_hat;
end


for i=1:k
Delsaboot(:,i)= bootstrp(nboot,@median,delsarmse(1:floor(size(x,1)/12),i));
end


figure
ylimm = [-0.05 0.55];
subplot(2,2,1)
boxplot(E(:,IND),'labels',labelsb(IND), 'labelorientation','inline')
h=findobj(gca,'tag','Outliers'); 
delete(h)
ylim(ylimm);
title('First Order Sobol \eta_i')

subplot(2,2,2)
boxplot(delta(:,IND),'labels',labelsb(IND), 'labelorientation','inline')
h=findobj(gca,'tag','Outliers'); 
delete(h) 
ylim(ylimm);
title('Borgonovo \delta_i')

subplot(2,2,3)
boxplot(T(:,IND),'labels',labelsb(IND), 'labelorientation','inline')
h=findobj(gca,'tag','Outliers'); 
delete(h) 
ylim(ylimm);
as=gca;
as.YTick=[0 0.1 0.2 0.3 0.4 0.5];
title('Kuiper index \beta^{KU}_i')

subplot(2,2,4)
boxplot(Delsaboot(:,IND),'labels',labelsb(IND), 'labelorientation','inline')
title('Median DELSA_i')
h=findobj(gca,'tag','Outliers'); 
delete(h) 
ylim(ylimm);
