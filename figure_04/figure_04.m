%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the convergence of the estimates for FUSE-016 model,
% increasing sample size from 200 to 9548
% Four sensitivity measures are considered:
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
addpath('./figure_04')

modelnr = 016;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);
delsarmse = csvread(['fuse',num2str(modelnr,'%03i'),'_delsa_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);

%% Convergence plot
rng(1234)

GDD = zeros(40,k);GDT = zeros(40,k);GDE = zeros(40,k);
N = [linspace(200,470,10), round(linspace(500,n,30))];

 for i = 1:length(N)
     M=30;
     nn = N(i);
    [b,d,t,e,w]=betaKS3(x(1:nn,:),y(1:nn),M); %B 2015 mehtod
    GDD(i,:) = d;
    GDT(i,:) = t;
    GDE(i,:) = e;
 end
 
ND = floor(N/(k+1));
GDC = zeros(40,k);
for ii=1:length(N)
    GDC(ii,:)=median(delsarmse(1:ND(ii),:));
end

%%
ylimm = [0,0.65];
refline1 = [600,600]; %eta
refline2 = [900,900]; %delta
refline3= [900,900]; %kuiper
refline4 = [600,600]; %delsa
figure
subplot(2,2,1)
hold on
a=1:11;
 for j = INDM
plot(N,GDE(:,j),'LineWidth',1.5)
 end
 
for j = INDR
    plot(N,GDE(:,j), 'LineStyle','--','color',[0.5 0.5 0.5])
end
hold off
line(refline1, [0 0.55], 'LineWidth',1,'LineStyle',':','color',[0.5 0.5 0.5]);
ylim(ylimm);xlabel('computational cost');
ylabel('\eta_i');title(['First Order Sobol'])

subplot(2,2,2)
hold on
a=1:11;
 for j =INDM
plot(N,GDD(:,j),'LineWidth',1.5)
 end
 
for j = INDR
    plot(N,GDD(:,j), 'LineStyle','--','color',[0.5 0.5 0.5])
end
hold off
line(refline2, [0 0.55], 'LineWidth',1,'LineStyle',':','color',[0.5 0.5 0.5]);
ylim(ylimm);xlabel('computational cost');
ylabel('\delta_i');title(['Borgonovo \delta'])

subplot(2,2,3)
hold on
a=1:11;
 for j = INDM
plot(N,GDT(:,j),'LineWidth',1.5)
 end
 
for j = INDR
    plot(N,GDT(:,j), 'LineStyle','--','color',[0.5 0.5 0.5])
end
hold off
line(refline3, [0 0.55], 'LineWidth',1,'LineStyle',':','color',[0.5 0.5 0.5]);
ylim(ylimm);xlabel('computational cost');
ylabel('\beta^{KU}_i');title(['Kuiper Index'])

subplot(2,2,4)
hold on
a=1:11;
 for j = INDM
plot(N,GDC(:,j),'LineWidth',1.5)
 end
 
for j = INDR
    plot(N,GDC(:,j), 'LineStyle','--','color',[0.5 0.5 0.5])
end
hold off
line(refline4, [0 0.55], 'LineWidth',1,'LineStyle',':','color',[0.5 0.5 0.5]);
ylim(ylimm);xlabel('computational cost');
ylabel('Median DELSA_i');title(['Median DELSA'])