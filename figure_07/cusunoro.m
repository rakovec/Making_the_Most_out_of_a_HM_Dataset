function [z,c]=cusunoro(x,y,param_names)
% CUSUNORO Shows CUSUNORO and 2nd moment CUSUNORO plots
%[z,v,c]=cusunoro(x,y,param_names)

% written by elmar.plischke@tu-clausthal.de
[n,k]=size(x);
[xr,indx]=sort(x);
yr=y(indx);
Ey=mean(y);
w=(y-Ey).^2;
Vy=sum(w);
c=(0:n)/n;
% CUSUNORO curve
z=[zeros(1,k);cumsum(yr-Ey)/ sqrt(n*Vy)];
if (nargin==3)
subplot(1,1,1);%subplot(1,2,1);
h=plot(c,z,'LineWidth',1.5);
if(nargin==3)
 legend(h,param_names{1:end-1});
 title(param_names{end});
else
    
legend({'TIMEDELAY','AXV\_BEXP','FRACTEN','MAXWATR\_1','MAXWATR\_2','PERCRTE','PERCEXP','BASERTE','QB\_POWR','LOGLAMB','TISHAPE'},'Location','EastOutside');
%legend( {'MAXWATR\_1','MAXWATR\_2','FRACTEN','PERCRTE','PERCEXP','BASERTE','QB\_POWR','AXV\_BEXP','LOGLAMB','TISHAPE','TIMEDELAY'},'Location','EastOutside');

%legend(cellstr(num2str((1:k)')),'Location','EastOutside')   
end
xlabel('cdf of inputs')
ylabel('cumulative sums of normalized output')

end
% 2nd moment CUSUNORO
% wr=w(indx);
% Ew=mean(w);
% Vw=sum((w-Ew).^2);
% v=[zeros(1,k);cumsum(wr-Ew)/sqrt(n*Vw)];
% subplot(1,2,2);
% plot(c,v);
% ylabel('cumulative sum for second moments');
% xlabel('cdf of inputs')
% if(nargin==3), title(param_names{end}); end
end