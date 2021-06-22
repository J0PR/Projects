%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Martyna Marczak
% Date: 31.12.2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab script generates cycle of the German real GDP, seasonally
% adjusted, in the framework of a structural time series model.
% It includes a function to perform boostrap in the way described by
% Stoffer and Wall (1998) and also a function to detect auxiliary
% residuals.
%
% The script uses the Toolbox SSMMatlab by Victor Gomez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

data = xlsread('data\PROJECTDATA.xls');
data(any(isnan(data)'),:) = [];

y = data(:,1);
ly = length(y);
Y=[];                    %matrix for regression variables
nreg=0;
npr=12;                  %number of forecasts
freq = 4;                % quarterly data
bg_year = 1970; bg_per = 1; 
datei = cal(bg_year,bg_per,freq);

% Specify components, initial values and fixed parameters

comp.level= [1  0 0];
comp.slope= [1  .005   NaN];
comp.irreg= [1  .1  NaN];
comp.cycle= [1 .1   NaN];
twopi=2*pi;
comp.cyclep=[0.9 twopi/40; NaN NaN];
comp.cycleb=[twopi/60. twopi/6.];   

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
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr
%Note that the standard errors are divided by the concentrated parameter
%(sqrt(result.sigma2c))

%create estimated model
[X,Z,G,W,T,H,ins,ii,ferror] = pr2usm(xvf,xf,str);


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
% tsplot(recrs, datei, 'GDP(STSM) St. Residuals');
% pause


%residual diagnostics
tname = 'GDP';
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
plotres([],[],[],[],[],1.96,'residuals',1,0,[],0,[],infr,1,1);

%file for output
fname='results\GDP_STSM.txt';
fid = fopen(fname,'w');
% fid=1;

%print estimation results
nreg=0;  
printusmer(fid,datei,tname,yor,y,ny,lam,str,result,nreg,nbeta);

%print residual diagnostics
printres(fid,infr);

%close external file
if fid ~= 1
 fclose(fid);
end

% Kalman smoothing
[KKP,PT,a,b]=scakfs(y,X,Z,G,W,T,H,ins,ii);

% Plot smoothed cycle
cyc = dispcomp(KKP,str,'cycle',datei,'GDP Smoothed Cycle');
%saveas(gcf, 'results\GDP_STSM','fig');

% % Create vector of times
% bg_times = (bg_year+(bg_per/freq));
% end_times = (bg_year+(bg_per/freq)) + (ly/freq) - (1/freq);
% times = (bg_times:(1/freq):end_times)';
% 
% % Write the cycle along with the times vector into the Excel file
% xlswrite('results\GDP_STSM.xls',[times cyc]);


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
fname='results\Aux_res_GDP.txt';
fid = fopen(fname,'w');
tname='GDP';

MDS1 = MDS1*conp;     

testauxres1(fid,tname,DS1,MDS1,str,result,dist,cw,1,2.3); 

%close external file
if fid ~= 1
 fclose(fid);
end

%close figures
closefig

%***********************************************************
%                  Bootstrap

nboot = 2; beta=result.h;        
% % the following is a logical array index to select the state vector 
% % components to be bootstrapped
% comp=logical([0 0 1 0]);     %cycle. 
% %If comp=logical([1 0 1 0]), then trend and cycle, etc.
comp='cycle';
Cycstar = bootssmcomp(y,X,Z,G,W,T,H,beta,str,nboot,comp);

