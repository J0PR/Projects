%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% US Industrial Production Index, 1946-I, 2011-III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

data=xlsread('c:\tst\series\PROJECT_US_MAN_RAW'); 

% GDP real wage: no intervention
% y = data(:,1);      %1946-I, 2011-III
y = data(57:end,1);      %1960-I, 2011-III
ly = length(y);
npr=12;                  %number of forecasts
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


%outliers detected by TRAMO; 1946-I, 2011-III
%   57 AO    ( 1 1960)
%  117 LS    ( 1 1975)
%   54 AO    ( 2 1959)
% Y=zeros(ly,3);           %matrix for regression variables
% nreg=3;
% Y(57,1)=1.;
% Y(117:end,2)=ones(ly-116,1);
% Y(54,3)=1.;
%outliers detected by TRAMO;    1960-I, 2011-III
%  61 LS    ( 1 1975)
Y=zeros(ly,2);           %matrix for regression variables
nreg=2;
Y(1,1)=1.;
Y(61:end,2)=ones(ly-60,1);

% Specify components, initial values and fixed parameters

comp.level= [1  0 0];
comp.slope= [1  0.005  NaN];
comp.seas=  [2  0.1 NaN];
comp.irreg= [1  .1  NaN];
comp.cycle= [1 0.1 NaN];
comp.conout='cycle';
twopi=2*pi;
comp.cyclep=[0.9 twopi/20; NaN NaN];
comp.cycleb=[twopi/40. twopi/6.];

comp.freq=freq;
comp.datei=datei;

if npr > 0
     mpr=npr; npr=0;
else
     mpr=0;
end

% Put univariate structural time series model into state space form
[str,ferror] = suusm(comp,y,Y,npr); 


if ferror > 0
return
end

% Estimate the state space model
[result,str] = usmestim(y,str);    

xvf=result.xvf;      % estimated parameters
xf=result.xf;        % fixed parameters

[X,Z,G,W,T,H,ins,ii,ferror] = pr2usm(xvf,xf,str);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Reduced form 
% %
% %estimation of the reduced form using the HR method
% x=[]; nfreq=1;
% yd=diffest(y,[],freq,0,1,1,0,0); %differenced series 
% [strv,ferror] = estvarmaxpqrPQR(yd,x,nfreq,[2 7 0],[0 0 0],0,1,1);
% %strv.phis3 contains the estimation of the cycle AR
% %strv.thetas3 contains the estimation of the MA part
% 
% %reduced form of the estimated structural model
% [nalpha,junk]=size(T); [ng,mg]=size(G);
% phiu(:,:,1)=eye(nalpha); phiu(:,:,2)=-T;
% thetau(:,:,2)=Z; np=nalpha+2;
% [phie,thetae,kro,ierror] = pright2leftcmfd(phiu,thetau,np);
% th=pmatmul(phie,G) + pmatmul(thetae,H); conp = result.sigma2c;
% cov=pmmulbf(th,th);   cov=cov.*conp;
% [nc,mc,pc]=size(cov); ninit=floor(pc/2)+1; Lp=cov(:,:,ninit:end);
% [Omega,Theta,ierror,iter,normdif]=pmspectfac(Lp,30,1e-4);
% %phie contains the AR part, Theta contains the MA part except Theta(0)=1.
% %
% % end of reduced form
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %compute recursive residuals
%beta=result.h;        
% %beta is the estimated regression vector. It will be considered fixed in 
% %the next function
%[Xtf,Ptf,recrs,recr,srecr,t1,A1,P1,KG]=scakfff(y,X,Z,G,W,T,H,ins,ii,beta);
% %recrs contains the standardized residuals
%t1 is the time at which the collapsed filter starts (no diffuse elements
%in the Kalman filter)
%A1 is x_{t1|t1-1}
%P1 is mse(A1)
% %KG contains the stack of the Kalman gains for the collapsed filter


% % Plot standardized recursive residuals (recrs)
% figure
% tsplot(recrs, datei, 'GWage(STSM) St. Residuals');
% pause


%residual diagnostics
tname = 'PR';
yor = y;
lam = 1;
e=result.e; F=result.F;
Ss=result.Ss; Ff=result.Ff; ne=length(e);    %residual sum of squares
Pevf=result.Pevf;                            %prediction error variance
% disp('standard error (finite sample)')
SPevf=result.SPevf;
ny=length(y); 
pvar=str.pvar; nr=length(pvar); 
X=str.X; [junk,nbeta]=size(X); 
W=str.W; [junk,nbetaw]=size(W); nbeta=nbeta+nbetaw;
ndrs=ne+nbeta;
lagl=min(36,max([floor(.2*ny) 3*freq 10]));
infr = rescomp(e,lagl,nr,Ss,Pevf,SPevf,Ff,ndrs,nbeta);  

%plot residual diagnostics  
% plotres([],[],[],[],[],1.96,'residuals',1,0,[],0,[],infr,1,1);

%file for output
fname='results\PR.txt';
fid = fopen(fname,'w');
% fid=1;
%print estimation results
% printusmerM(fid,datei,tname,yor,y,ny,lam,str,result,nreg,nbeta);
printusmer(fid,datei,tname,yor,y,ny,lam,str,result,nreg,nbeta); 

%print residual diagnostics
printres(fid,infr);

%close external file
if fid ~= 1
 fclose(fid);
end

% 
% Kalman smoothing
[KKP,PT,a,b]=scakfs(y,X,Z,G,W,T,H,ins,ii);


Cc = dispcomp(KKP,str,'cycle',datei,'PR Smoothed Cycle');
cyc = Cc(:,1);

return

% % Plot smoothed cycle
% cyc = KKP(:,6);
% figure
% tsplot(cyc, datei, 'PR Smoothed Cycle');
% pause
% saveas(gcf, 'results\PR','fig');

%save cycle 
save results\cpr.dat cyc -ascii -double


% % Create vector of times
% bg_times = (bg_year+(bg_per/freq));
% end_times = (bg_year+(bg_per/freq)) + (ly/freq) - (1/freq);
% times = (bg_times:(1/freq):end_times)';
% 
% % Write the cycle along with the times vector into the Excel file
% xlswrite('results\GWage_STSM.xls',[times cyc]);


%***********************************************************
%                  Bootstrap
% 
nboot = 2;
comp = 'cycle';
%comp=logical([0 0 1 0]);     %cycle.
%If comp=logical([1 0 1 0]), then trend and cycle, etc.
beta=result.h; 
Cycstar = bootssmcomp(y,X,Z,G,W,T,H,beta,str,nboot,comp); 
%*************************************************************************

%*************************************************************************
%                   Tests with auxiliary residuals
%
 conp = result.sigma2c;   % concentrated out variance
%
% Computation with the function smoothgen.m
% 
% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = D_t*epsilon_t
% Hence, U_t = C_t = 0

% nalpha = ii(1);
[nalpha,junk]=size(T); 
mucd = 1+nalpha;
U = [];
% C = zeros(ly*mucd,nalpha);
% D = repmat([G;H],ly,1); 
C = zeros(mucd,nalpha);
D = [G;H];

[DS1,MDS1,hd1,Md1] = smoothgen(y,X,Z,G,W,T,H,ins,ii,mucd,U,C,D);

dist = {'irreg'};
cw= 1.96;
fname='results\Aux_res_PR.txt';
fid = fopen(fname,'w');
tname='PR';

MDS1 = MDS1*conp;     

testauxres1(fid,tname,DS1,MDS1,str,result,dist,cw,1,3.0); 
%close external file
if fid ~= 1
 fclose(fid);
end

%close figures
closefig

