function [INDM,INDR,IND,labels,labelsb]=getparaindex(modelnr,k)

fid = fopen(['fuse',num2str(modelnr,'%03i'),'_parameters_base.csv']);
labels = textscan(fid,'%s',1);
fclose(fid);

labels = cellstr(labels{1});
labels =strsplit(labels{1},',');
labelsb = strrep(labels,'"','');
labels =  strrep(labelsb,'_','\_');

if modelnr == 016
    INDM=[11 8 3];
    INDR =[1 2 4 5 6 7 9 10];
    IND = [INDM,INDR];
elseif modelnr == 014
    INDM=[13 10 6];
    INDR =setdiff([1:k],INDM);
    IND = [INDM,INDR];
elseif modelnr == 072
    INDM=[11 8 5 3];
    INDR =setdiff([1:k],INDM);
    IND = [INDM,INDR];
    elseif modelnr == 160
    INDM=[12 11 5 3];
    INDR =setdiff([1:k],INDM);
    IND = [INDM,INDR];
    elseif modelnr == 170
    INDM=[14 11 7];
    INDR =setdiff([1:k],INDM);
    IND = [INDM,INDR];
end

end