%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots D-scatter plot for FUSE-160 model
%
% written by Xuefei LU
% 2016 Dec
% email: xuefei.lu@unibocconi.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% import data
currentpath = cd('..');
addpath('./data_input')
addpath('./figure_S3')

modelnr = 160;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
sensitivity_x= csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_sensitivities_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(sensitivity_x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);

%% D-Scatter plot
figure
vnames = labels(IND);
for i=1: k
xx = x(:,IND(i));
yy=sensitivity_x(:,IND(i));
ind = (yy <=0);
if (strcmp(labels(IND(i)),  'TIMEDELAY')) ind = (yy <=0.1);end
subplot(ceil(k/3+.5),3,i)
hold on
scatter(xx(ind),yy(ind),3,[0    0.4470    0.7410],'filled')
scatter(xx(~ind),yy(~ind),3,[0.8500    0.3250    0.0980],'filled')
if (strcmp(labels(IND(i)),  'TIMEDELAY')) ylim([-15,2]);end
if (strcmp(labels(IND(i)),   'AXV\_BEXP')) ylim([-10,10]);end
if (strcmp(labels(IND(i)),  'PERCRTE'))  ylim([-0.05 0.05]); end
if (strcmp(labels(IND(i)),  'BASERTE'))  ylim([-0.006 0.002]); end
if (strcmp(labels(IND(i)),  'LOGLAMB'))  ylim([-0.5 1]); end
if (strcmp(labels(IND(i)),   'TISHAPE'))  ylim([-1 0.5]); end
if (strcmp(labels(IND(i)),   'QB\_POWR'))  ylim([-2 1]); end
 xlim([min(xx), max(xx)]);
xlabel([labels(IND(i))]);
ylabel('Partial derivative')
hold off
end
 