clearvars
uqlab
rng(1234) % Set the random seed

currentpath = cd('..');
addpath('./data_input')
addpath('./figure_08')

modelnr = 016;
x = csvread(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv'],1);
y = csvread(['fuse',num2str(modelnr,'%03i'),'_raw_rmse.csv'],1);

[n,k] = size(x);
[INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k);
vnames = labels(IND);
x=x(:,IND);
%% PCE METAMODEL
nk=2000; %using N=2000 points to train PCE subroutine
XT=x(1:nk,:);
YT=y(1:nk);

myInput = CreateInput(vnames);
clear metaopts
metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.ExpDesign.X = XT;
metaopts.ExpDesign.Y = YT;
metaopts.Degree = 1:10; 
metaopts.Input = myInput;
myPCE =  uq_createModel(metaopts);

%% Calculate second order Sobol
rng(1234567)
PSobol.Type = 'Sensitivity';
PSobol.Method = 'Sobol';
PSobol.Sobol.Order = 2;
PSobol.SampleSize=100000;
PSobolAnalysis = uq_createAnalysis(PSobol);
FirstOrderSobol = PSobolAnalysis.Results.FirstOrder;
TotalOrderSobol = PSobolAnalysis.Results.Total;
%uq_display(PSobolAnalysis);

resultind = PSobolAnalysis.Results.VarIdx{2};
SIOv=PSobolAnalysis.Results.AllOrders{2};
[Sv,Sind] = sort(SIOv,'descend');%2nd Sobol estimates
Pind =resultind(Sind(1:5),:); % top 5 index
PCEs = Sv(1:5);

%% save data
save('PCEresult.mat','Pind', 'PCEs','FirstOrderSobol','TotalOrderSobol' )