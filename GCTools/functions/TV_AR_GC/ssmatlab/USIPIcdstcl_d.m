%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% US Industrial Production Index, 1946-I, 2011-III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

data=xlsread('c:\tst\series\PROJECT_US_MAN_RAW'); 

% y = data(:,1);      %1946-I, 2011-III
y = data(57:end,1);      %1960-I, 2011-III
ly = length(y);
freq = 4;                % quarterly data
bg_year = 1960; bg_per = 1; 
datei = cal(bg_year,bg_per,freq);

lag=24; cw=1.96;  
tname={'US Industrial Production Index'};
for dr=0:1     
 for ds=0:1
  c0=sacspacdif(y,tname,dr,ds,freq,lag,cw);
  pause
 end
end
closefig


% Model identified by TRAMO is (1,1,0)(0,1,1)_4 without mean
% outliers detected by TRAMO;    1960-I, 2011-III
%  61 LS    ( 1 1975)
Y=zeros(ly,1);           %matrix for regression variables
nreg=1;
Y(61:end,1)=ones(ly-60,1);


%estimate model using HR method 
%model is
dr=1; p=1; q=0;
ds=1; P=0; Q=1;
yd=diffest(y,[],freq,0,dr,ds,0,0); %difference series 
Yd=diffest(Y,[],freq,0,dr,ds,0,0);
%preliminary estimation of beta using OLS
beta=mulols(yd,Yd); ydc=yd-Yd*beta;
%Hannan-Rissanen method to estimate ARIMA parameters
x=[];
[strv,ferror] = estvarmaxpqrPQR(ydc,x,freq,[p dr q],[P ds Q],0,1,1);  

%setup model
phi(:,:,1)=1; Phi(:,:,1)=1;  th(:,:,1)=1;  Th(:,:,1)=1; 
phi(:,:,2)=strv.phis3(:,:,2); 
Th(:,:,2)=strv.thetas3(:,:,freq+1);
Sigma=strv.sigmar3; 

%create structure and put model into state space form
[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq);   

%estimate model using exact likelihood method
result=varmapqPQestim(yd,str,Yd);
%estimated and fixed parameters
xvf=result.xvf; xf=result.xf;  
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr
%residual diagnostics
e=result.e;      %white noise residuals
ff=result.ff;    %vector of nonlinear functions
nbeta=0; %length of regression vector
ndrs=length(yd); 
Ss=e'*e; Ff=ff'*ff;    %residual sum of squares
conp=result.sigma2c; sconp=sqrt(conp);
lagl=4*freq;
infr = rescomp(e,lagl,length(xvf),Ss,conp,sconp,Ff,ndrs,nbeta);
%plot residual diagnostics  
plotres([],[],[],[],[],1.96,'residuals',1,0,[],0,[],infr,1,1);
closefig
%print residual diagnostics
%file for output
fname='results\PRbp.txt';
fid = fopen(fname,'w');
% fid=1;
printres(fid,infr);
%close external file
if fid ~= 1
 fclose(fid);
end

%*************************************************************************
%                CANONICAL DECOMPOSITON
%  
% set up estimated model
[phif,thf,Phif,Thf,Lf,ferror] = pr2varmapqPQ(xvf,xf,str); 

%phi = [-0.500677731003569   1.]; Th  = [-0.855929516763095   1.];
% stand. dev. of resid (sconp) = 0.0167476578257468

 %update structure str with the estimated parameters
 %Note that the residual covariance matrix is divided by the concentrated 
 %parameter (result.sigma2c = conp).
 Sigmaf=(Lf*Lf').*conp; 
 [str,ferror] = suvarmapqPQ(phif,thf,Phif,Thf,Sigmaf,freq); 


% set up trend-cycle and seasonal polynomials for the canonical
% decomposition
s=freq;   
[phir,phis,thr,ths,phirst] = arima2rspol(phif,Phif,thf,Thf,s,dr,ds);   

% perform canonical decomposition
% [thrc,sigma2r,thsc,sigma2s,thtc,sigma2t,sigma2i,ierrcandec] = ...
%  candec(phir,phis,thr,ths,phirst,dr,ds,sconp); 
[compcd,ierrcandec] = candec(phir,phis,thr,ths,phirst,s,dr,ds,sconp);
sigma2i=compcd.itvar;
if sigma2i < 0
 disp('irregular spectrum negative')
%  return
 disp('model is changed. sigma2i is made equal to zero.')
 disp('this is a provisional solution to the negative ')
 disp('irregular spectrum problem')
 sigma2i=0.;
end
%   
% This canonical decomposition gives rise to the following
%
%                   Models for the components  
%
% thrc =  0.1523   -0.9679   -0.1202    1.0000
% phir = -0.5007    2.0014   -2.5007    1.0000
% sigma2r = 0.314624552730715
% thsc =  0.0200300618841744  1.02902270418928   1.46898735203652   1.
% phis =  1.                  1.                 1.                 1.
% sigma2s = 0.00197541841872233
% sigma2i = 0.0956411020166517
% stand. dev. of resid (sconp) = 0.0167476578257468


% Models for the components obtained by SEATS (Canonical decomposition),
% assuming that the model has been estimated by TRAMO
% % Model estimated by TRAMO
% % phi=[-.49555 1]; Th=[-.86020 1]; stand. dev. of resid.: 0.1661D-01
%
%                   MODELS FOR THE COMPONENTS
%
%  TREND-CYCLE NUMERATOR
%    1.00000000000000      -0.120915955672721      -0.968885573465440     
%   0.152030382207280     
%  TREND-CYCLE DENOMINATOR
%      1.0000    -2.4956     1.9911    -0.4956
%  INNOV. VAR. (*)     0.31544
% 
%  SEAS. NUMERATOR
%    1.00000000000000        1.46717885612683        1.02591246661574     
%   1.787979416833018E-002
%  SEAS. DENOMINATOR
%      1.0000     1.0000     1.0000     1.0000
%  INNOV. VAR. (*)     0.00187
% 
%  IRREGULAR
%  VAR.     0.09674
%
% (*)   IN UNITS OF VAR(A)
% STANDARD DEVI. OF RESID=  0.1661D-01
%


% put the canonical decomposition model y_t = Y_t*beta + p_t + s_t + r_t + 
% i_t into state space form (Akaike state space form for each component)
% create structure
npr=0;
[X,Z,G,W,T,H,ins,ii,strc,ferror] = sucdm(compcd,y,Y,str,npr);     

% Kalman smoothing for the canonical model
[KKP,PT,a,b]=scakfs(y,X,Z,G,W,T,H,ins,ii);




% stochastic trend
% trend = KKP(:,1);
Cc = dispcomp(KKP,strc,'trendcycle');  
trend = Cc(:,1);
vnames=strvcat('PR','PR trend');  
figure
tsplot([y trend],datei,vnames);
pause
% closefig
% %save trend 
% save results\tpr.dat trend -ascii -double

%%
% ***********************************************************************
%                        LOW-PASS FILTERS
%
% design of low-pass filter to obtain a smooth trend and a cycle. The
% filter is specified giving Lambda and Di. The sine But. filter is the
% Hodrick-Prescott filter. The filter will be applied to the canonical
% trend.
Lambda=1600; Di=2;
% select from the two following lines the tangent or the sine But. filter
% by uncommenting or commenting the appropriate lines.
% [compf,ferror]=dtanbut([],[],[],Di,[],Lambda); 
[compbst,ferror]=dsinbut([],[],[],Di,[],Lambda); 
% plot gain function of both low-pass filters
figure
ggsintanbut([],[],[],compbst.Di,compbst.Thetac)
pause
% put the model y_t = Y_t*beta + sp_t + c_t + s_t + r_t + i_t into Akaike 
% state space form, where sp_t is the smooth trend and c_t is the cycle, 
% both obtained from the previous p_t by application of the low pass 
% filter. 
[X,Z,G,W,T,H,ins,ii,strc,ferror] = sucdmpbst(compcd,compbst,y,Y,str,npr);    
nsp=length(compbst.den) + length(compcd.ptden)-1;

% Kalman smoothing for the canonical plus bsttrend model
[KKP,PT,a,b]=scakfs(y,X,Z,G,W,T,H,ins,ii);


% smooth trend and cycle obtained by applying the low-pass filter. Plot the
% cycle first
Cc = dispcomp(KKP,strc,{'trend','cycle'},datei,2,'PR bstcycle');
bsttrend = Cc(:,1);
bstcycle = Cc(:,2); 
% bsttrend = KKP(:,1);
% bstcycle = KKP(:,nsp+1);
% vnames='PR bstcycle';  
% figure
% tsplot(bstcycle,datei,vnames);
% pause


% plot both the smooth and the canonical trend
vnames=strvcat('PR trend','PR bsttrend');  
figure
tsplot([trend bsttrend],datei,vnames);
pause
% closefig


%************************************************************************
%                    BAND-PASS FILTER 
%
% design of band-pass filter to obtain a well defined cycle and a 
% relatively smooth trend. Frequencies are expressed divided by pi. The
% filter will be applied to the canonical trend. 

D(1)=.1; D(2)=.1; xp1=.0625; xp2=.3; xs=.4;


% Tangent band-pass filter
[compbp,ferror]=dbptanbut(D,xp1,xp2,xs);
% 
% plot gain function of the tangent band-pass filter
figure
ggbptanbut(D,xp1,xp2,xs,compbp.Di,compbp.Alph,compbp.Lambda)
pause
closefig



% put the model y_t = Y_t*beta + sp_t + c_t + s_t + r_t + i_t into Akaike 
% state space form, where sp_t is the smooth trend and c_t is the cycle, 
% both obtained from the previous p_t by application of the band-pass 
% filter. 
% 

[X,Z,G,W,T,H,ins,ii,strc,ferror] = sucdmpbp(compcd,compbp,y,Y,str,npr); 
nsp=length(compbp.den) + length(compcd.ptden)-1;

% Kalman smoothing for the canonical plus bpcycle model
[KKP,PT,a,b]=scakfs(y,X,Z,G,W,T,H,ins,ii);  


% cycle obtained by applying the band-pass filter and smooth trend obtained
% as difference between the canonical trend and the cycle. Plot the cycle
% first.
Cc = dispcomp(KKP,strc,{'trend','cycle'},datei,2,'PR bpcycle');
bptrend = Cc(:,1);
bpcycle = Cc(:,2);
% bptrend = KKP(:,1);
% bpcycle = KKP(:,nsp+1);
% vnames='PR bpcycle';  
% figure
% tsplot(bpcycle,datei,vnames);
% pause


% plot both the smooth and the canonical trend
vnames=strvcat('PR trend','PR bptrend');  
figure
tsplot([trend bptrend],datei,vnames);
pause

% compare band-pass cycle with cycle obtained by applying the low-pass 
% filter
figure
vnames=strvcat('PRbtcycle','PR bpcycle');  
tsplot([bstcycle bpcycle],datei,vnames);
pause

closefig

% % save band-pass cycle
% save results\PRUSbp.dat bpcycle -ascii -double


