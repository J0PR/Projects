%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 8.2 of Reinsel (1997), pp. 292-298
%
% Series are: 1) the in-phase current, 2) the out-of-phase current and 3)
% the frequency of the voltage generated, expressed as deviations from
% nominal values. The number of observations is n=100. The data have been
% multiplied by 10. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y=load('c:\ssm_matlab\data\power-turbo.dat'); 

lag=20; cw=1.96; freq=1; ds=0; 
tname={'In-phase current deviations','Out-of-phase current deviations',...
       'Frequency of voltage deviations'};
for i=1:3      
 for dr=0:1
  c0=sacspacdif(y(:,i),tname(i),dr,ds,freq,lag,cw);
  pause
 end
end
closefig


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


%Although lagsopt=5, we estimate a VAR of order 4 as in Reinsel p.295
test=1; lags=4;
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

%To estimate the conditional ML estimate of Sigma, we estimate the model
%again using the Hannan-Rissanen method, that in this case coincides with
%the conditional method. 
seas=1; x=[]; hr3=1;
[strv,ferror] = estvarmaxpqrPQR(y,x,seas,[4 0 0],[0 0 0],hr3); 

disp('Covariance matrix of residuals estimated by conditional ML:')
disp(strv.sigmar2) 


%fix unsignificant paramaters to zero and estimate again
strv.phi(1:2,3,2)=zeros(2,1); strv.phi(1:2,3,3)=zeros(2,1);
strv.phi(1:2,3,4)=zeros(2,1); strv.phi(1:2,3,5)=zeros(2,1);
strv.nparm=strv.nparm-8;
strv=mhanris(y,x,seas,strv,hr3);

disp('Covariance matrix of residuals estimated by conditional ML')
disp('for the restriced model:')
disp(strv.sigmar2)
pause

lag=6; ic=1; recrs=strv.resid2;
str=mautcov(recrs,lag,ic); 

disp('Lag-1 residual correlation matrix:')
disp(str.r(:,:,1))
pause

%estimate a VARMA(4,1) by the Hannan-Rissanen method 
clear strv
hr3=0; finv2=1;
[strv,ferror] = estvarmaxpqrPQR(y,x,seas,[4 1 0],[0 0 0],hr3,finv2);

%impose the restriction that all of the (1,3) and (2,3) elements are zero  
%and estimate again
strv.phi(1:2,3,2)=zeros(2,1); strv.phi(1:2,3,3)=zeros(2,1);
strv.phi(1:2,3,4)=zeros(2,1); strv.phi(1:2,3,5)=zeros(2,1);
strv.theta(1:2,3,2)=zeros(2,1);
strv.nparm=strv.nparm-10;
strv=mhanris(y,x,seas,strv,hr3,finv2);

%estimate the model using the conditional method
[xvfc,strc,ferrorc]=mconestim(y,x,strv);  
pause



%impose additional constraints as in Reinsel (1997)
strv.phi(2,1,2)=0; strv.phi(3,2,2)=0; strv.phi(2,1,3)=0; strv.phi(3,2,3)=0;
strv.phi(3,3,3)=0;
strv.phi(3,3,4)=0; strv.phi(3,2,5)=0;
strv.theta(3,:,2)=zeros(1,3); strv.theta(2,1,2)=0; strv.theta(1,2,2)=0;
strv.nparm=strv.nparm-12;
strv=mhanris(y,x,seas,strv,hr3,finv2);

%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(y,x,strv);
pause

%estimate again to improve a little over the results
[xvfc,strc,ferrorc]=mconestim(y,x,strc);
pause


%estimate using exact ML
%setup model
Phi=eye(3); Th=eye(3);
phi=strc.phiscon; th=strc.thetascon(:,:,1:2); Sigma=strc.sigmarcon;  

%create regression variable for the mean
Y=eye(3); 
%fix insignificant parameters to zero and set up model  
phi(1:2,3,2)=zeros(2,1); phi(1:2,3,3)=zeros(2,1);
phi(1:2,3,4)=zeros(2,1); phi(1:2,3,5)=zeros(2,1);
theta(1:2,3,2)=zeros(2,1);
phi(2,1,2)=0; phi(3,2,2)=0; phi(2,1,3)=0; phi(3,2,3)=0;
phi(3,3,3)=0;
phi(3,3,4)=0; phi(3,2,5)=0;
th(3,:,2)=zeros(1,3); th(2,1,2)=0; th(1,2,2)=0;
freq=1;
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq);  
str.phin(1:2,3,2)=zeros(2,1); str.phin(1:2,3,3)=zeros(2,1);
str.phin(1:2,3,4)=zeros(2,1); str.phin(1:2,3,5)=zeros(2,1);
str.thn(1:2,3,2)=zeros(2,1);
str.phin(2,1,2)=0; str.phin(3,2,2)=0; str.phin(2,1,3)=0; str.phin(3,2,3)=0;
str.phin(3,3,3)=0;
str.phin(3,3,4)=0; str.phin(3,2,5)=0;
str.thn(3,:,2)=zeros(1,3); str.thn(2,1,2)=0; str.thn(1,2,2)=0;
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
lag=6; ic=1;
str=mautcov(recrs,lag,ic);  
disp('sample autocorrelations signs:')
disp(str.sgn)
pause

%estimate a VARMAX(4,0,4) by the Hannan-Rissanen method. Var y_t3 is the
%output and y_1t and y_2t are the inputs. See p. 297 in Reinsel (1997)
x=y(:,1:2); yo=y(:,3);
[strv,ferror] = estvarmaxpqrPQR(yo,x,seas,[4 0 4],[0 0 0],hr3,finv2);
strv.phi(1,1,3)=0; strv.phi(1,1,4)=0; 
strv.gamma(1,1,1)=0; strv.gamma(1,2,1)=0;
strv.gamma(1,2,2)=0; strv.gamma(1,2,3)=0; strv.gamma(1,2,5)=0;
strv.nparm=strv.nparm-7;
strv=mhanris(yo,x,seas,strv,hr3,finv2);
disp('estimated phi parameters:')
disp(strv.phis3)
disp('estimated gamma parameters:')
disp(strv.gammas3)
pause

%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(yo,x,strv);
disp('estimated phi parameters:')
disp(strc.phiscon)
disp('estimated gamma parameters:')
disp(strc.gammascon)
pause

% Y=1;
% [xvf,strx,ferror]=mexactestim(yo,x,strc,Y);
