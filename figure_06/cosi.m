function Si=cosi(x,y,M,gfx,vnames)
%% COSI Calculation of sensitivity indices from given data.
%     SI = COSI(X,Y) returns the sensitivity indices for input arguments X 
%     and output arguments Y (per line) based upon a discrete cosine 
%     transformation.
%%
%     Written by Elmar Plischke, elmar.plischke@tu-clausthal.de


 [n,k]=size(x);
 [nn,kk]=size(y);
 if nn~=n, error('Input/output sizes mismatch!'), end

 [xr,index]=sort(x);
 if kk==1
% sort output
    yr=y(index);
 else
    yr=zeros(n,k*kk);
    for i=1:kk
        z=y(:,i);
        yr(:,(i-1)*k+(1:k))=z(index);
    end
 end
 
 %% frequency selection
if (nargin==2) || (isempty(M))
  M=max(ceil(sqrt(n)),3);
 fprintf('COSI: Using %d coefficients.\n',M);
end 
% consider M terms
d=zeros(1,n);
d(1+(1:M))=1;

%% Compute transformation
allcoeff=dct(yr);

% transformation is orthogonal, so by Parseval's Theorem
V = sum(allcoeff(2:end,:).^2);
Vi= sum(allcoeff(1+(1:M),:).^2);
Si= Vi./V;
%disp('mono?'), -allcoeff(2,:)
%% 
if nargin==5
%vnames = {'MAXWATR\_1','MAXWATR\_2','FRACTEN','PERCRTE','PERCEXP','BASERTE','QB\_POWR','AXV\_BEXP','LOGLAMB','TISHAPE','TIMEDELAY'};
%vnames={'TIMEDELAY','AXV\_BEXP','FRACTEN','MAXWATR\_1','MAXWATR\_2','PERCRTE','PERCEXP','BASERTE','QB\_POWR','LOGLAMB','TISHAPE'};
 for i=1:k
  if(k>1), subplot(ceil(k/3+.5),3,i); end
  yhat=zeros(n,1);
  yhat(1:(M+1))=allcoeff(1:(M+1),i);
  plot(x(index(:,i),i),y(index(:,i)),'.',x(index(:,i),i),idct(yhat),'-','LineWidth',2);
  xlim([min(x(index(:,i),i)), max(x(index(:,i),i))]);
  xlabel([vnames(i)]);%xlabel(['Input {',num2str(i),'}']);
  ylabel('raw\_rmse');
  %title(gfx);
 end
end
%%%
 if kk>1, Si=reshape(Si',k,kk)'; end
 return
