 function [d,s,p]=deltafast2(x,y,M,gfx,vnames)
% DELTAFAST a quick method of computing delta moment independent measure
% DELTA=DELTAFAST(X,Y)
% [DELTA,THETA,LAMBDA]=DELTAFAST(X,Y) Kullback-Leibler/Shannon
% Information/L2 pdf density
%delta/ K-L/  L2 pdf density
%%
[n,k]=size(x);
% Epanechnikov with variance 1
Kernel=@(x)3/(4*sqrt(5))*max(1-(x.^2/5),0);
% Box kernel
%Kernel=@(x)(1-abs(x/sqrt(3))>0)/(2*sqrt(3));

% Number of partitions
if(nargin<3)||isempty(M), M=24; end
% Numerical noise cutoff (simple Kolmogorov-Smirnov)
Cutoff = 0.7; % 1.31; % 1.5; % 0..2 % Missing Factor 2 ???

%% output stats
%ys=sort(y);
miny=min(y);
maxy=max(y);
L=miny-.04*(maxy-miny);
U=maxy+.04*(maxy-miny);

%% transform to unbounded support
% no trafo
Stretch=@(y,l,u)y;
Squeeze=@(z,l,u)z;
% log trafo
%Stretch=@(y,l,u)log(y-l)-log(u-y);
%Squeeze=@(z,l,u)(exp(z)*u+l)./(exp(z)+1);
% atanh trafo
%Stretch=@(y,l,u)atanh((2*y-(l+u))/(u-l));
%Squeeze=@(z,l,u)(tanh(z)*(u-l)+(l+u))/2;
% probit trafo
%Stretch=@(y,l,u)norminv((y-l)/(u-l));
%Squeeze=@(z,l,u)l+normcdf(z)*(u-l);

ty=Stretch(y,L,U); 

% Transform to Gaussian
%[ty,iy]=sort(y);yr(iy)=1:n;
%ty(iy)=-sqrt(2)*erfinv(1-(2*(1:n)-1)'/n);
%% work with transformed data
medy=median(ty);
iqry=median(abs(medy-ty)); % interquartile range estimator
% bandwidth estimate (rule of thumb)
stdy=min(std(ty),iqry/0.675);
h=1.2*stdy*((4/(3*n))^(1/5)); % 1.2 compensates subset selection
% construct interpolation points
z1=linspace(min(ty)-2*h,medy-iqry, 25);
z2=linspace(medy-iqry,medy+iqry,52);
z3=linspace(medy+iqry,max(ty)+2*h,25);
z=[z1,z2(2:end-1),z3];
%z=linspace(L,U,110);
l=length(z);
% back-trafo interpolation points for gfx
tz=Squeeze(z,L,U);
%% kernel density matrix
%W=Kernel( (repmat(z,n,1)-repmat(ty,1,l))/h)/h;
W=Kernel( bsxfun(@minus,z,ty)/h)/h;
%% unconditional density 
densy=mean(W);

%% conditional densities for partitioned data
if(nargin==5),cols=jet(M);plotCols=ceil(sqrt(k));clf;end
Sm=zeros(k,M);
Tm=zeros(k,M);
Lm=zeros(k,M);
nm=zeros(k,M);

%% keep only W from the partition
[xr,indxx]=sort(x);
for i=1:k
   xr(indxx(:,i),i)=1:n; % ranks (no ties)
end
%%
for j=1:M
   indx= ((j-1)*n/M <xr) & (xr <= j*n/M);
   nm(:,j)=sum(indx); % no ties: always same nr. of realizations
   for i=1:k
   % conditional density
   densc=mean(W(indx(:,i),:));
   % L1 separation of densities (using Scheffe Thm)
   Sm(i,j)=trapz(z,max(densy-densc,0)); % all entries <1
   % Kullback Leibler
   tt=densc.*(log(densc)-log(densy));tt(densc==0)=0;
   Tm(i,j)=trapz(z,tt);
   % Power score - L2 pdf density
   Lm(i,j)=trapz(z,(densy-densc).^2);
   if(nargin==5)
    % vnames={'TIMEDELAY','AXV\_BEXP','FRACTEN','MAXWATR\_1','PERCEXP','MAXWATR\_2','PERCRTE','BASERTE','QB\_POWR','LOGLAMB','TISHAPE'};
     subplot(ceil(k/plotCols)+1,plotCols,i);
     if(j==1)||(j==M),plot(tz,densy,'k','LineWidth',3);hold on;end
     plot(tz,densc,'Color',[0 0.447 0.741])%cols(j,:));
      title([vnames(i)]);
      xlabel('rmse');
      ylim([0,1])
   end
   end
% Clear noise
Sm(Sm<Cutoff.*sqrt(1/n+1./nm))=0;
if(nargin==5) 
 subplot(ceil(k/plotCols)+1,1,ceil(k/plotCols)+1);
 %plot(1:M,Sm.*nm/n)
 plot(cumsum(nm(1,:))/n-1/(2*M),Sm);%L1 separation of densities (using Scheffe Thm)
 hold on
 %plot(cumsum(nm(1,:))/n-1/(2*M),Tm,':'); %Kullback Leibler
 %plot(cumsum(nm(1,:))/n-1/(2*M),2*sqrt(Lm),'--');
 title('Separation of Conditional Densities, L^1 norm')
 xlabel('Empirical cdf of parameters');
 hold off
 legend([vnames(1:5)],'Location','EastOutside');
 a=axis;a(1)=0;a(2)=1;axis(a);
end 
 d=sum(Sm.*nm,2)'/n;
 s=sum(Tm.*nm,2)'/n;
 p=sum(Lm.*nm,2)'/n;
end
