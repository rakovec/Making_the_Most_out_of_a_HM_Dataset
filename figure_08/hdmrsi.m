function [Sij,R2]=hdmrsi(u,y,M,gfx)
% HDMRSI estimates first and second order effects.
% It uses a polynomial regression with shifted Legendre polynomials 
% as suggested by Rabitz et al. for Random-Sampling High-Dimensional
% Model Representation
%  SIJ=HDMRSI(U,Y) 
%  [SIJ,R2]=HDMRSI(U,Y)
%  SIJ=HDMR(U,Y,[MAX1,MAX2]) additionally specifies
%  the maximal polynomial order for additive and quadratic models
%  SIJ=HDMR(U,Y,[],'Graphics Title') illustrates the model fits
% U is assumed uniformly [0,1]

%Sij=hdmrsi(uniformer(x),y,[6,4])
%matrix Sji with joint effects in the upper, interaction effects in the lower off-diagnal part
% degree 6 for first order effects, 4 for 2nd order effects


% Written by elmar.plischke@tu-clausthal.de
%%
[n,k]=size(u);
[nn,kk]=size(y);
if nn~=n, error('Input/output sizes mismatch!'), end
if kk~=1, error('Only single output supported.'); end

%%
if(nargin<3 || isempty(M))
    M=[4,4];
end
%% Lanczos: App.Ana., Table IV
Pn = [-1    2      0      0       0       0       0        0       0       0 0;
       1   -6      6      0       0       0       0        0       0       0 0;
      -1   12    -30     20       0       0       0        0       0       0 0;
       1  -20     90   -140      70       0       0        0       0       0 0;
      -1   30   -210    560    -630     252       0        0       0       0 0;
       1  -42    420  -1680    3150   -2772     924        0       0       0 0;
      -1   56   -756   4200  -11550   16632  -12012     3432       0       0 0;
       1  -72   1260  -9240   34650  -72072   84084   -51480   12870       0 0;
      -1   90  -1980  18480  -90090  252252 -420420   411840 -218790   48620 0;
       1 -110   2970 -34320  210210 -756756 1681680 -2333760 1969110 -923780 184756];
%%
Sij=zeros(k,k);
EY=mean(y);
VY=sum((y-EY).^2);
%% design matrix: constant factor
D1=ones(n,1); 
%% design matrix: additive factors
%Di=[sqrt(3)*(2*u-1), sqrt(5)*(6*u.^2-6*u+1),sqrt(7)*(20*u.^3-30*u.^2+12*u-1)];
Di=zeros(n,k*M(1));
%% for visualisation
pts1d=30;pts2d=10;
v=linspace(0,1,pts1d);Fi=zeros(pts1d,M(1)); 
[wx,wy]=meshgrid(linspace(0,1,pts2d));ww=[wx(:),wy(:)];
Gi=zeros(pts2d^2,M(1)*2);
% pre-compute the phi-terms
for i=1:M(1)
% too slow
%    % normalised shifted Legendre polynomials, cf. Lanczos: App.Ana.
%%    Di(:,(i-1)*k+(1:k))= hypergeom(sqrt(2*i+1)*[-i,i+1],1,1-u);
% using table lookup
    phi=Pn(i,(i+1):-1:1)*sqrt(2*i+1);
% using symbolic calculus
% Matlab 7.10 only likes i<=3
%    phi=sqrt(2*i+1)*sym2poly(hypergeom([-i,i+1],1,'1-x'));
    Di(:,(i-1)*k+(1:k))= polyval(phi,u);
    Fi(:,i)            = polyval(phi,v);
    Gi(:,2*i-[1,0])    = polyval(phi,ww);
end
%% fit one-parameter model
for(i=1:k)
    DD=Di(:, i+(0:(M(1)-1))*k);
    D=[D1,DD];
    beta=D\y; 
    %yhat=D*beta; % for first order effects
    if(nargin==4)
        figure(1);
        subplot(floor(k/2+.5),2,i);       
        % plot(u(:,i),y,'*',u(:,i),yhat,'.',v,[ones(30,1),Fi]*beta,'k-o');
        plot(u(:,i),y,'*');hold on;
        plot(v,[ones(pts1d,1),Fi]*beta,'-','Color',[0,.5,0],'LineWidth',2);hold off;
        ylabel('Output');
        xlabel(['Input parameter u_{' num2str(i) '}']);
    end
    Sij(i,i)=n*sum( beta(2:end).^2)/VY; %sum( (yhat -EY).^2)/VY;
end 
%% design matrix: interactions
 % after some quiet hours blind-foldedly single-stepping through gui_hdmr ...
 z=ones(M(1),1)*(1:M(1));
 z=1i*z'+z;
phi_combinations=[real(z(:)),imag(z(:))];
%% for 2D gfx
Gij=[ones(pts2d^2,1), Gi(:,1:2:(2*M(1)-1)), Gi(:,2:2:(2*M(1)))];
for c=phi_combinations'
    if(sum(c)<=M(2))
        Gij = [ Gij, Gi(:, 2*c(1)-1).*Gi(:, 2*c(2)) ];
    end
end
%%
Dij=[];
for(i=1:k)
    for (j=(i+1):k)
        Dj=[];
        for c=phi_combinations'
            if(sum(c)<=M(2))
                Dj = [ Dj, Di(:, i+(c(1)-1)*k).*Di(:, j+(c(2)-1)*k) ];
            end
        end
        Dij=[Dij, Dj];
        D=[D1,Di(:, i+(0:(M(1)-1))*k), Di(:, j+(0:(M(1)-1))*k),Dj];
        beta=D\y;
        yhat=D*beta;
        Sij(i,j)=sum( (yhat -EY).^2)/VY;
        Sij(j,i)=n*sum(beta( (2+2*M(1)):end).^2)/VY;
        if nargin==4
            figure(2)
            surf(wx,wy,reshape(Gij*beta,pts2d,pts2d));
            shading interp
            hold on
            plot3(u(:,i),u(:,j),y,'*'); %,u(:,i),u(:,j),yhat,'+');
            xlabel(['u_{' num2str(i) '}']);
            ylabel(['u_{' num2str(j) '}']);
            zlabel('Output')
            hold off;
            pause
        end
    end
end
%%
R2=sum(sum(tril(Sij)));
%%
if(0)
%% alternative: fit additive model
D=[D1,Di];
beta=D\y;
yhat=D*beta;
R2=sum( (yhat -EY).^2)/VY % goodness of fit
%% now activate only relevant coefficients for parameter i
for i=1:k
    betai=beta;
    indx=1:M(1)*k;
    indx(i+(0:(M(1)-1))*k)=[];
    betai( 1+indx)=0;
    yhati=D*betai;
    Si(i)=sum( (yhati -mean(yhati)).^2)/VY;
%    Sj(i)=sum( (yhati -EY).^2)/VY;
%    Sk(i)=sum(betai(2:end).^2)*n/VY;
end
Si %,Sj,Sk
%%
end
%%
end
