%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the Second order interatcion effects for FUSE-016 model
%The Second Order Sobol indices are estiamted using 3 subroutines: 
%HDMR : cheap version of Ziehn Tomlinson, by Elmar Plischke
%LASI: by Elmar Plischke
%PCE: by UQlab
%
% written by Xuefei LU
% 2016 Dec
% email: xuefei.lu@unibocconi.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% import data
currentpath = cd('..');
addpath('./data_input')
addpath('./figure_08')

modelnr = 016;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);

%% important parameter indices and labels
[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);
vnames = labels(IND);
x=x(:,IND);

%% HDMR
transU=@(X)bsxfun(@rdivide,bsxfun(@minus,X,min(X)),max(X)-min(X));
u=transU(x); %transform inputs into uniform
[Sij,R2]=hdmrsi(u,y,[6,4]); 
% Find the 5 largest 2nd order effect
Sint=tril(Sij) ;
Sint(logical(eye(size(Sint)))) = 0;
[sortedValues,sortIndex] = sort(Sint(:),'descend');    
maxIndex = sortIndex(1:5);  %# Get a linear index into A of the 5 largest values
[Hind(:,1),Hind(:,2)]=ind2sub(size(Sij),maxIndex);
HDMRs = sortedValues(1:5);
                                     
%% lasi
opts=struct('M',7); 
[SiL,SijL]=lasi(u,y,opts);
Z=SijL; 
for i=1:size(x,2), E=zeros(size(SijL));E(i,:)=1;E(:,i)=1;Z=Z-E*SiL(1,i); end
SZ=tril(Z); %joint effect, =triu(Sij)
Sint=tril(SZ); 
Sint(logical(eye(size(Sint)))) = 0;
[sortedValues,sortIndex] = sort(Sint(:),'descend');    
maxIndex = sortIndex(1:5) ; 
[Lind(:,1),Lind(:,2)]=ind2sub(size(Sij),maxIndex);
LASIs= sortedValues(1:5);

%% Import PCE result
load('PCEresult.mat')


%% figure_08a
figure
subplot(3,1,1)
for i=1:5
vnamesK(i)=strcat(vnames(min(Pind(i,:))),'\newline',vnames(max(Pind(i,:))));
end
bar(PCEs,0.9,'FaceColor',[0.9290 0.6940 0.1250])
set(gca,'YLim',[0,0.12],'fontsize',9, 'XTick', 1:5, 'XTickLabel', vnamesK);
title(['PCE£¬N=2000'])
ylabel('S_{i,j}')
ylim([0,1])
subplot(3,1,2)
for i=1:5
vnamesH(i)=strcat(vnames(min(Hind(i,:))),'\newline',vnames(max(Hind(i,:))));
end
bar(HDMRs,0.9,'FaceColor',[0 0.4470 0.7410])
set(gca,'YLim',[0,0.12],'fontsize',9,'XTick', 1:5, 'XTickLabel', vnamesH); 
title(['HDMR£¬N=',num2str(n)])
ylabel('S_{i,j}')
ylim([0,1])
subplot(3,1,3)
for i=1:5
vnamesL(i)=strcat(vnames(min(Lind(i,:))),'\newline',vnames(max(Lind(i,:))));
end
bar(LASIs,0.9,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'YLim',[0,0.12],'fontsize',9, 'XTick', 1:5, 'XTickLabel', vnamesL);
title(['LASI£¬N=',num2str(n)])
ylabel('S_{i,j}')
ylim([0,1])

%% figure_08b
figure
s=80;
c=[get(gca,'ColorOrder');
    0.4 0.8 1;
    1 1 0.2;
    1 0 0.6;
    0.2 1 0.6;
    0 0.4 0.6;
    0 0.6 0.8;
    0.4 0 0];    
sp = ['o','d','s','^','v','>','.','.','.','.','+', '*', 'x','.'] ;
hold on
for i=1:6
scatter(FirstOrderSobol(i),TotalOrderSobol (i), s,c(i,:),'filled',sp(i))
end
for i=7:k
scatter(FirstOrderSobol(i),TotalOrderSobol (i,:), s,c(i,:),sp(i))
end
hold off
ml=min(get(gca,'XLim'),get(gca,'YLim'));
line(ml,ml,'linewidth',1,'color',[0.5,0.5,0.5],'LineStyle','--');
legend(vnames,'Location','EastOutside')
xlabel('First Order Sobol')
ylabel('Total Order Sobol')
