%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots cusunoro plot for FUSE-160 model
%
% written by Xuefei LU
% 2016 Dec
% email: xuefei.lu@unibocconi.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% import data
currentpath = cd('..');
addpath('./data_input')
addpath('./figure_S2')

modelnr = 160;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);

%% cusunoro plot
figure
[z,c]=cusunoro(x(:,IND),y);
bb=round(linspace(1,n,35));
LS=char([ {'--o','--^','--s','--','--x'},repelem({':'}, 1, k-5)]);
hold on
for ii=1:k
plot(c(bb),z(bb,ii),LS(ii,:),'LineWidth',1); %
end
hold off
xlabel('cdf of inputs')
ylabel('cumulative sums of normalized output')
legend(labels(IND),'Location','EastOutside');