function [myInput] = CreateInput(vnames)
%%  create input object for PCE 

k=length(vnames);

for i = 1 : k
    Input.Marginals(i).Type = 'Uniform' ;
    if (strcmp(vnames(i),  'TIMEDELAY'))
        Input.Marginals(i).Parameters =[0.1,2];end
    if (strcmp(vnames(i),   'AXV\_BEXP'))
        Input.Marginals(i).Parameters =[0.001, 3];end
    if (strcmp(vnames(i),  'PERCEXP'))
        Input.Marginals(i).Parameters =[1, 20]; end
    if (strcmp(vnames(i),  'MAXWATR\_1'))
        Input.Marginals(i).Parameters =[50, 500]; end
    if (strcmp(vnames(i),  'MAXWATR\_2'))
        Input.Marginals(i).Parameters =[25, 250]; end
    if (strcmp(vnames(i),  'FRACTEN'))
        Input.Marginals(i).Parameters =[0.05,0.95]; end
    if (strcmp(vnames(i),  'FRCHZNE'))
        Input.Marginals(i).Parameters =[0.05,0.95]; end
    if (strcmp(vnames(i),  'PERCRTE'))
        Input.Marginals(i).Parameters =[0.01, 1000]; end
    if (strcmp(vnames(i),  'FRACLOWZ'))
        Input.Marginals(i).Parameters =[0.05, 0.95]; end
    if (strcmp(vnames(i),  'BASERTE'))
        Input.Marginals(i).Parameters =[0.001, 1000]; end
    if (strcmp(vnames(i),   'QB\_POWR'))
        Input.Marginals(i).Parameters =[1, 10]; end
    if (strcmp(vnames(i),  'LOGLAMB'))
        Input.Marginals(i).Parameters =[5, 10]; end
    if (strcmp(vnames(i),   'TISHAPE'))
        Input.Marginals(i).Parameters =[2, 5]; end
    %160
    if (strcmp(vnames(i),   'RTFRAC1'))
        Input.Marginals(i).Parameters =[0.05, 0.95]; end
    if (strcmp(vnames(i),  'FPRIMQB'))
        Input.Marginals(i).Parameters =[0.05, 0.95]; end
    if (strcmp(vnames(i),  'PERCFRAC'))
        Input.Marginals(i).Parameters =[0.05, 0.95]; end
    if (strcmp(vnames(i),  'QBRATE\_2A'))
        Input.Marginals(i).Parameters =[0.01, 0.25]; end
    if (strcmp(vnames(i),  'QBRATE\_2B'))
        Input.Marginals(i).Parameters =[0.01, 0.25]; end
    
end

myInput = uq_createInput(Input);