function [Si,Sij]=lasi(x,y,opts_in)
%% LASI Calculation of sensitivity indices from given data.
%     SI = LASI(X,Y) returns the sensitivity indices for input arguments X 
%     and output arguments Y (per line) based upon a harmonic regression.
%     [SI,SIJ]=LASI(X,Y) additionally returns the paired group effects
%     SI = LASI(X,Y,OPTS) computes the indices using the options structure
%           OPTS containing the fields
%                   'M' [6]             - max. higher harmonic
%                   'BackTrafo' [@(x)x] - Input trafo to unit cube, 
%                                          use a string for empirical cdf
%                   'Gfx' ['Test']      - title for optional graphics
%%
%     Written by Elmar Plischke, elmar.plischke@tu-clausthal.de

opts=struct('M',6,'BackTrafo',@(x)x,'Gfx','Test','PlotCols',2,'Gfx3D','');
if nargin>2 && isstruct(opts_in)
        members=fieldnames(opts);
        for i=1:length(members)
            o=members{i};
            if isfield(opts_in,o), opts.(o)=opts_in.(o);end
        end
end        
[n,k]=size(x);
[nn,kk]=size(y);
if nn~=n, error('Input/output sizes mismatch!'), end
%%
EY=mean(y);
VY=var(y);
Si=zeros(2,k);
%Sj=zeros(1,k);
% BackTrafo is either empirical cdf or simulated cdf
if(isa(opts.BackTrafo,'function_handle'))
 u=opts.BackTrafo(x);
else
 % do empirical cdf
 u=zeros(n,k);
 for i=1:k
     [xs,js]=sort(x(:,i));
     u(js,i)=((1:n)-.5)/n; % (0:(n-1))/(n-1);
 end
end
 if(any(min(u)<0) && any(max(u)>1)), 
     error('LASI: Transformed values outside [0,1]'); 
 end
 D=ones(n,opts.M+1)/sqrt(n);
 if(nargout>1)
  Ds=zeros(n,opts.M+1,k);
 end
%%
for (i=1:k)
%% Construct design matrix for harmonic (cosine) fit
    for (j=1:opts.M)
     D(:,j+1)=sqrt(2/n)*cos( pi*u(:,i)*j);
    end
%    if(isa(opts.Group,'numeric'))
%        if(i==opts.Group(1)) D1=D; end
%        if(i==opts.Group(2)) D2=D; end
%    end
   if(nargout>1)
       Ds(:,:,i)=D;
   end
%  Least squares regression    
    c=D\y;
%  First order sensitivity indicator    
    Si(1,i)=sumsqr(c(2:end))/(n-1)/VY;
% First order sensitivity - alternative formulation
    yhat=D*c;
    Si(2,i)=sumsqr(yhat-EY)/(n-1)/VY;
% Graphics output
    if(~isempty(opts.Gfx))
%         if(k>1), subplot(floor(k/opts.PlotCols+.5),opts.PlotCols,i); end
%         [xs,js]=sort(x(:,i));
%         plot(u(:,i),y,'.',u(js,i),yhat(js),'k-','LineWidth',2,'MarkerSize',1);
%         xlabel(['u_' num2str(i)]);ylabel('y');
%         title(opts.Gfx)
    end
%%    
end
%% Test for second order
    if nargout>1
    if(~isempty(opts.Gfx3D))
        figure
        l=ceil(sqrt(k*(k-1)/2));
        c=1;
    end
    Sij=zeros(k,k);
%    if(isa(opts.Group,'numeric'))
%        Sij=zeros(2,1);
    for i=1:(k-1);
      for j=(i+1):k
        DD=xmul(Ds(:,:,i),Ds(:,:,j));
        cc=DD\y;
        % sum of squares coefficients
        Sij(i,j) = sumsqr(cc(2:end))/(n-1)^2/VY;
        yhat=DD*cc;
        if(~isempty(opts.Gfx3D))
            subplot(l,l,c);c=c+1;
          plot3(u(:,i),u(:,j),y,'.',u(:,i),u(:,j),yhat,'k*','MarkerSize',3);
          xlabel(['u_{' num2str(i) '}']);
          ylabel(['u_{' num2str(j) '}']);
        end
        % sum of squares model predictions
        Sij(j,i) = sumsqr(yhat-EY)/(n-1)/VY;
    end,end,end
end

% row-wise kronecker product
function C=xmul(A,B)
[n,k]=size(A);
[m,l]=size(B);
if(n~=m), error('Size mismatch.'); end
C=zeros(n,k*l);
for (i=1:n)
    C(i,:)=kron(B(i,:),A(i,:));
end
end
