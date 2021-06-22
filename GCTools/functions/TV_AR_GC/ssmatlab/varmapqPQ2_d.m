%Example of estimation of a VARMA(p,q)(P,Q)_s model in which we fix some
%parameters. The initial estimates are obtained using the Hannan-Rissanen
%method. For the estimation, the conditional method and the exact methods
%are used. 
%Series is simulated series from Tiao and Box (1981)
% -- THIS IS A SIMULATED VECTOR ARMA(1,1) EXAMPLE
% --INPUT VARI ARE CONST, PHI, THETA, MEAN, SIGMA. NCOLS 1,2,2,1,2.
% -- 1.20  0.7 0.0  0.5 0.4  0.0  1.0 1.2
% -- 2.20  0.4 0.9  0.0 0.5  0.0  1.2 4.0
% --END OF DATA
% --MTSMODEL NAME IS ARMA11. SERIES ARE S1,S2. MODEL IS  @
% --  (1 - PHI*B)SERIES=CONST + (1 - THETA*B)NOISE.  
%

clear
y=load('data\vf_classarma.dat'); x=[];

% define model. 
freq=1; 
npr=12;                  %number of forecasts
%copy npr in mpr and make npr zero for estimation
if npr > 0
 mpr=npr; npr=0;  
else
 mpr=0;
end

lag=20; cw=1.96; ds=0; dr=0;
tname={'TB simulated ARMA(1,1) 1','TB simulated ARMA(1,1) 2'};
for i=1:2      
  c0=sacspacdif(y(:,i),tname(i),dr,ds,freq,lag,cw);
  pause
end
closefig


%estimate model using HR method 
seas=1;
[strv,ferror] = estvarmaxpqrPQR(y,x,seas,[1 1 0],[0 0 0],0,1,1);   

%estimate using the conditional method
[xvfc,strv,ferrorc]=mconestim(y,x,strv); 


%setup model
Phi=eye(2); Th=eye(2);
phi=strv.phiscon; th=strv.thetascon; Sigma=strv.sigmarcon;  

%create structure and put model into state space form
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq);
%create regression variable for the mean
Y=eye(2); 

%estimate model using the exact method
result=varmapqPQestim(y,str,Y);   

%check t-values
disp('t-values: ')
disp(result.tv)
disp('press any key to continue')
pause


%fix two insignificant parameters to zero and set up model again
phi(1,2,2)=0.; th(2,1,2)=0.;   
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq); 
str.phin(1,2,2)=0; str.thn(2,1,2)=0; 
[str,ferror] = fixvarmapqPQ(str); 
%reestimate model
result=varmapqPQestim(y,str,Y);   

%estimated and fixed parameters
xvf=result.xvf; xf=result.xf;  
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr

%create estimated model
[phif,thf,Phif,Thf,Lf,ferror] = pr2varmapqPQ(xvf,xf,str);   

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
%see +-. signs
disp(str.sgn)
pause
disp(str.sgnt)
pause

%compute forecasts
if mpr > 0
 %hb, Mb, A and P are in structure result. Here, hb is the vector of
 %regression estimates and Mb is the matrix of standard errors. A is the
 %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the 
 %end of the sample and P is the matrix of standard errors. 
 hb=result.h; Mb=result.H; A=result.A; P=result.P;

 npr=mpr;
 Xp=Y;  
 Wp=[];
 cw=1.96;
 m=2;                           %number of series
 [pry,mypr,alpr,malpr]=ssmpred(npr,m,A,P,Xp,Z,G,Wp,T,H,hb,Mb);
 conp=result.sigma2c;  
 spry=zeros(m,npr); sconp=sqrt(result.sigma2c);
 for i=1:npr
  spry(:,i)=sqrt(diag(mypr(:,:,i)))*sconp;
 end
 opry=pry; ospry=spry;
 %plot forecasts 
 tname='var1';
 out.pry=pry(1,:); out.spry=spry(1,:); 
 out.opry=opry(1,:); out.ospry=ospry(1,:); out.y=y(:,1); 
 out.yor=y(:,1); out.ny=length(y(:,1)); out.npr=npr; out.cw=cw; 
 out.tname=tname;
 lam=1;                   %lam=0, logs are taken; =1, no logs are taken
 out.lam=lam; out.s=freq;
 pfctsusm(out);
 tname='var2';
 out.pry=pry(2,:); out.spry=spry(2,:); 
 out.opry=opry(2,:); out.ospry=ospry(2,:); out.y=y(:,2); 
 out.yor=y(:,2); out.ny=length(y(:,2)); out.npr=npr; out.cw=cw; 
 out.tname=tname;
 lam=1;                   %lam=0, logs are taken; =1, no logs are taken
 out.lam=lam; out.s=freq;
 pfctsusm(out);
end