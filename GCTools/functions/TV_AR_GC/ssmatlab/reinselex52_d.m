%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 5.2 of Reinsel (1997), pp. 170-174
%
% Series are: 1) weekly production schedule figures and 2) billing figures.
% The number of observations is n=100.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y=load('c:\ssm_matlab\data\Weeklyproshed.dat'); 

[ny,s]=size(y);

lag=20; cw=1.96; freq=1; ds=0; 
tname={'weekly production schedule figures','billing figures'};
for i=1:2      
 for dr=0:1
  c0=sacspacdif(y(:,i),tname(i),dr,ds,freq,lag,cw);
  pause
 end
end
closefig

%sample correlation matrices 
lag=10; ic=1;  
str=mautcov(y,lag,ic);
disp('Sample Correlation matrices (lags 1-5)')
str.r(:,:,1:5)
pause
disp('Sample Correlation matrices (lags 6-10)')
str.r(:,:,6:10)
pause


%Preliminary VAR analysis
prt=1;  
minlag=0; maxlag=6;
lagsoptlr = lratiocr(y,maxlag,minlag,prt);
disp('lagsoptlr')
disp(lagsoptlr)
pause

crt='bic';
lagsoptbic = infcr(y,maxlag,minlag,crt,prt);  
disp('lagsoptbic')
disp(lagsoptbic)
pause

crt='aic';
lagsoptaic = infcr(y,maxlag,minlag,crt,prt);  
disp('lagsoptaic')
disp(lagsoptaic)
pause

%Estimate a VAR of order 5  
test=1; lags=5;
res = var_est(y,lags,test);

disp('estimated VAR coefficient matrices:')
disp(res.betavar')

disp('t-values:')
disp(res.tvvar')

disp('Estimated covariance matrix of residuals:')
disp(res.sigmar)

disp('Granger causality prob.:')
disp(res.fprob)
pause

lag=12; ic=1; nr=s^2*lags;
str=mautcov(res.resid,lag,ic,nr);

disp('sample autocorrelations signs:')
disp(str.sgnt)
pause

disp('p-values of Q statistics:')
disp(str.pval)
pause


%estimate the model with two variables  
%First, estimate a VARMAX(4,1,0) model by the Hannan-Rissanen method. 
x=[]; hr3=0; finv2=1;  
[strv,ferror] = estvarmaxpqrPQR(y,x,freq,[4 1 0],[0 0 0],hr3,finv2); 

%Then, impose restrictions and estimate using the HR method again
strv.phi(1,1,2:5)=zeros(1,4);  strv.phi(1,2,2:5)=zeros(1,4);
strv.phi(2,1,2:3)=zeros(1,2); strv.phi(2,2,4:5)=zeros(1,2);
strv.theta(1,2,2)=0; strv.theta(2,1,2)=0;
strv.nparm=strv.nparm-14;
strv=mhanris(y,x,freq,strv,hr3,finv2);   
disp('estimated phi parameters:')
disp(strv.phis3)
pause
disp('estimated theta parameters:')
disp(strv.thetas3)
pause


%estimate using exact ML
%setup model
Phi=eye(2); Th=eye(2);
phi=strv.phis3; th=strv.thetas3(:,:,1:2); Sigma=strv.sigmar3;  

%create regression variable for the mean
Y=eye(2); 
%fix insignificant parameters to zero and set up model  
Sigma(2,1)=0; Sigma(1,2)=0;
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq); 
str.phin(1,1,2:5)=zeros(1,4);  str.phin(1,2,2:5)=zeros(1,4);
str.phin(2,1,2:3)=zeros(1,2); str.phin(2,2,4:5)=zeros(1,2);
str.thn(1,2,2)=0; str.thn(2,1,2)=0; 
str.Lhn(2)=0;
[str,ferror] = fixvarmapqPQ(str);  

%estimate model using the exact method
result=varmapqPQestim(y,str,Y);  

%estimated and fixed parameters
xvf=result.xvf; xf=result.xf;  
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr

%create estimated model
[phif,thf,Phif,Thf,Lf,ferror] = pr2varmapqPQ(xvf,xf,str);
Sigmar=result.Sigmar;
%t-values
tvf=result.tv; 
[phitvf,thtvf,Phitvf,Thtvf,Ltvf,ferror] = pr2varmapqPQ(tvf,xf,str);  

disp('estimated phi:')
disp(phif)
pause
disp('estimated th:')
disp(thf)
pause

disp('t-values of phi:')
disp(phitvf)
pause
disp('t-values of th:')
disp(thtvf)
pause

disp('Estimated covariance matrix of residuals:')
disp(result.Sigmar)

%compute recursive residuals
%Note that the residual covariance matrix is divided by the concentrated 
%parameter (result.sigma2c).
Sigmaf=Lf*Lf'; 
[strf,ferror] = suvarmapqPQ(phif,thf,Phif,Thf,Sigmaf,freq); 
%set up regression matrices
X=Y; W=[]; 
%set up system matrices
T=strf.T; Z=strf.Z; G=strf.G; H=strf.H;
%set up initial state
% [Sigma,ferror] = mlyapunov(T,H*H'); 
Sigma = dlyapsq(T,H);  Sigma=Sigma'*Sigma;
ins=Sigma; 
[nalpha,mf]=size(strf.T);
i=[nalpha 0 0 0]; s=str.s;

[Xt,Pt,g,M,initf,recrs]=scakff(y,X,Z,G,W,T,H,ins,i);
%plot recursive residuals
plot(recrs(:,1)), legend('recrs(:,1)'),pause
plot(recrs(:,2)), legend('recrs(:,2)'),pause
closefig
%compute autocovariance and autocorrelation matrices of rec. residuals
lag=12; ic=1; nr=length(xvf)-s*(s+1)/2+1;
str=mautcov(recrs,lag,ic,nr);  
disp('sample autocorrelations signs:')
disp(str.sgn)
pause

disp('p-values of Q statistics:')
disp(str.pval)
pause




%alignment
y1a=y(:,1); y1a(3:end)=y(1:98,1);
ya=[y1a(3:end) y(3:end,2)];

%Preliminary VAR analysis
prt=1;  
minlag=0; maxlag=6;
lagsoptlr = lratiocr(ya,maxlag,minlag,prt);
disp('lagsoptlr')
disp(lagsoptlr)
pause

crt='bic';
lagsoptbic = infcr(ya,maxlag,minlag,crt,prt);  
disp('lagsoptbic')
disp(lagsoptbic)
pause

crt='aic';
lagsoptaic = infcr(ya,maxlag,minlag,crt,prt);  
disp('lagsoptaic')
disp(lagsoptaic)
pause

%Estimate a VAR of order 3 as in Reinsel p.173
test=1; lags=3;
res = var_est(ya,lags,test);

disp('estimated VAR coefficient matrices:')
disp(res.betavar')

disp('t-values:')
disp(res.tvvar')

disp('Estimated covariance matrix of residuals:')
disp(res.sigmar)

disp('Granger causality prob.:')
disp(res.fprob)
pause

lag=12; ic=1;  nr=s^2*lags;
str=mautcov(res.resid,lag,ic,nr);

disp('sample autocorrelations signs:')
disp(str.sgnt)
pause

disp('p-values of Q statistics:')
disp(str.pval)
pause


%estimate a VARMA(2,1) by the Hannan-Rissanen method 
hr3=0; finv2=1; x=[];
[strv,ferror] = estvarmaxpqrPQR(ya,x,freq,[2 1 0],[0 0 0],hr3,finv2);

%estimate the model using the conditional method
[xvfc,strc,ferrorc]=mconestim(ya,x,strv);  
pause

disp('estimated phi:')
disp(strc.phiscon)
pause
disp('estimated th:')
disp(strc.thetascon)
pause

disp('t-values of phi:')
disp(strc.phitvcon)
pause
disp('t-values of th:')
disp(strc.thetatvcon)
pause

%fix all insignifican parameters to zero and estimate again
strv.phi(1,1:2,2)=zeros(1,2); strv.phi(1,1:2,3)=zeros(1,2);
strv.theta(2,1,2)=0;
strv.nparm=strv.nparm-5;
strv=mhanris(ya,x,freq,strv,hr3,finv2);

%estimate the model using the conditional method
[xvfc,strc,ferrorc]=mconestim(ya,x,strv);  
pause

disp('estimated phi:')
disp(strc.phiscon)
pause
disp('estimated th:')
disp(strc.thetascon)
pause

disp('t-values of phi:')
disp(strc.phitvcon)
pause
disp('t-values of th:')
disp(strc.thetatvcon)
pause

disp('estimated covariance matrix of residuals')
disp(strc.sigmarcon)
pause

lag=12; ic=1; nr=length(xvf)-s*(s+1)/2+1;
str=mautcov(strc.residcon,lag,ic,nr);

disp('sample autocorrelations signs:')
disp(str.sgnt)

disp('p-values of Q statistics:')
disp(str.pval)
pause

