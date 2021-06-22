%Example of estimation of an ARMA(p,q)(P,Q)_s model
%Series is airline series from Box and Jenkins (1976)
%

clear
% load data. 
y=load('data\bjsgairl.dat'); x=[];
yl=log(y);                      %transform series


lag=36; cw=1.96; freq=12;  
tname={'BJ Airline Passengers'};
for dr=0:1     
 for ds=0:1
  c0=sacspacdif(yl,tname,dr,ds,freq,lag,cw);
  pause
 end
end
closefig

yd=diffest(yl,[],12,0,1,1,0,0); %difference series 

% define model. Model is (0,0,1)(0,0,1)_12 for the differenced logged
% series
phi(:,:,1)=1; Phi(:,:,1)=1; 
th(:,:,1)=1;  Th(:,:,1)=1; 
%no mean in the model
Y=[];
npr=12;                  %number of forecasts
%copy npr in mpr and make npr zero for estimation
if npr > 0
 mpr=npr; npr=0;  
else
 mpr=0;
end

%estimate model using HR method 
[strv,ferror] = estvarmaxpqrPQR(yd,x,freq,[0 1 0],[0 1 0],0,1,1);  

%setup model
th(:,:,2)=strv.thetas3(:,:,2); Th(:,:,2)=strv.thetas3(:,:,freq+1);
Sigma=strv.sigmar3;


%create structure and put model into state space form
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq);   

%estimate model
result=varmapqPQestim(yd,str,Y);    

%estimated and fixed parameters
xvf=result.xvf; xf=result.xf;  
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr

%create estimated model
[phif,thf,Phif,Thf,Lf,ferror] = pr2varmapqPQ(xvf,xf,str);

%residual diagnostics
e=result.e;      %white noise residuals
ff=result.ff;    %vector of nonlinear functions
nbeta=0; %length of regression vector
ndrs=length(yd); 
Ss=e'*e; Ff=ff'*ff;    %residual sum of squares
conp=result.sigma2c; sconp=sqrt(conp);
lagl=3*freq;
infr = rescomp(e,lagl,length(xvf),Ss,conp,sconp,Ff,ndrs,nbeta);

%plot residual diagnostics  
plotres([],[],[],[],[],1.96,'residuals',1,0,[],0,[],infr,1,1);
closefig

%print residual diagnostics
%file for output
fname='results\bjsgairl.txt';
fid = fopen(fname,'w');
% fid=1;
printres(fid,infr);
%close external file
if fid ~= 1
 fclose(fid);
end

%compute forecasts of the logged differenced series
if mpr > 0
 %hb, Mb, A and P are in structure result. Here, hb is the vector of
 %regression estimates and Mb is the matrix of standard errors. A is the
 %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the 
 %end of the sample and P is the matrix of standard errors. 
 hb=result.h; Mb=result.H; A=result.A; P=result.P;

 npr=mpr;
 %set up system matrices for the estimated model
 %Note that the residual covariance matrix is divided by the concentrated 
 %parameter (result.sigma2c).
 Sigmaf=Lf*Lf'; 
 [strf,ferror] = suvarmapqPQ(phif,thf,Phif,Thf,Sigmaf,freq); 
 Z=strf.Z; G=strf.G; T=strf.T; H=strf.H;
 Xp=Y;  
 Wp=[];
 cw=1.96;
 m=1;                           %number of series
 [pry,mypr,alpr,malpr]=ssmpred(npr,m,A,P,Xp,Z,G,Wp,T,H,hb,Mb);
 spry=zeros(m,npr); sconp=sqrt(result.sigma2c);
 for i=1:npr
  spry(:,i)=sqrt(diag(mypr(:,:,i)))*sconp;
 end
 %obtain forecasts in the original scale using the log-normal
 %distribution
 opry=pry; ospry=spry;

 %plot forecasts 
 tname='bjsgairl (differenced and in logs)';
 out.pry=pry; out.spry=spry; out.opry=opry; out.ospry=ospry; out.y=yd; 
 out.yor=yd; out.ny=length(yd); out.npr=npr; out.cw=cw; out.tname=tname;
 lam=1;                   %lam=0, logs are taken; =1, no logs are taken
                          %in this case, since we work with the logged 
                          %series, lam=1.
 out.lam=lam; out.s=freq;
 pfctsusm(out);
 
end

%compute forecasts of the original series 
if mpr > 0
 npr=mpr;
 %set up system matrices for the estimated ARIMA model
 %Note that the residual covariance matrix is divided by the concentrated 
 %parameter (result.sigma2c).
 Sigmaf=Lf*Lf'; 
 %Differencing polynomial
 phifo(:,:,1)=1.; phifo(:,:,2)=-1.; Phifo(:,:,1)=1.; Phifo(:,:,2)=-1.;
 %MA polynomial
 thfo=thf; Thfo=Thf;
 [strfo,ferror] = suvarmapqPQ(phifo,thfo,Phifo,Thfo,Sigmaf,freq);
 %ARIMA model in state space form
 Z=strfo.Z; G=strfo.G; T=strfo.T; H=strfo.H; [ndelta,junk]=size(T);
 X=[]; W=[];
 %initial conditions for the Kalman filter
 i=[0 0 0 ndelta]; ins=eye(ndelta); chb=1;
 % [ins,i,ferror]=incossm(T,H,ndelta);
 %the whole initial vector is nonstationary. Otherwise use the previous
 %commented line
 
 %run Kalman filter
 [e,f,hb,Mb,A,P,qyy,R]=scakfle2(yl,X,Z,G,W,T,H,ins,i,chb);
 %hb is the vector of regression estimates and Mb is the matrix of standard 
 %errors. A is the estimated state vector, x_{t|t-1}, obtained with the 
 %Kalman filter at the end of the sample and P is the matrix of standard 
 %errors. 
 
 %forecasts
 [pry,mypr,alpr,malpr]=ssmpred(npr,m,A,P,Xp,Z,G,Wp,T,H,hb,Mb);
  spry=zeros(m,npr); sconp=sqrt(result.sigma2c);
 for i=1:npr
  spry(:,i)=sqrt(diag(mypr(:,:,i)))*sconp;
 end
 %obtain forecasts in the original scale using the log-normal
 %distribution
 lam=0;
 opry=pry; ospry=spry; 
 if lam == 0
  for i=1:npr
   opry(i)=exp(pry(i)+(spry(i)^2)/double(2.));
   ospry(i)=exp(double(2.)*pry(i)+spry(i)^2)*(exp(spry(i)^2)-double(1.));
  end
 end
 %plot forecasts 
 tname='bjsgairl';
 out.pry=pry; out.spry=spry; out.opry=opry; out.ospry=ospry; out.y=yl; 
 out.yor=y; out.ny=length(yl); out.npr=npr; out.cw=cw; out.tname=tname;
 out.lam=lam; out.s=freq;
 pfctsusm(out);
end

