%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Martyna Marczak
% Date: 31.12.2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab script generates a cycle for the German consumer real wage 
% series, seasonally adjusted. It uses the framework of structural time 
% series models. The series is assumed to have a slope intervention in the 
% first quarter of 2003 (an impulse).
% It includes a function to perform boostrap in the way described by
% Stoffer and Wall (1998).
%
% The script uses the Toolbox SSMMatlab by Victor Gomez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

data = xlsread('data\PROJECTDATA.xls');
data(any(isnan(data)'),:) = [];

y = data(:,2);
ly = length(y);

freq = 4;                % quarterly data
bg_year = 1970;
bg_per = 1; 
datei = cal(bg_year,bg_per,freq);

% Determine the observation number at the time of the intervention
i_year = 2003;           % starting year of the intervention
i_per = 1;               % starting period of the intervention
% observation number at the time of the intervention
in = ((i_year + i_per) - (bg_year + bg_per))*freq + 1;

% Incorporate the structural break into the equation for trend slope
% Intervention variable is in this case a pulse variable, i.e. it 
% takes the value 1 at the time point of the intervention and 0 otherwise.
%
% The state space model is given by:
% alpha_{t+1} = W_t * beta + T_t * alpha_t + H_t * eps_t
%         Y_t = X_t * beta + Z_t * alpha_t + G_t * eps_t
%
% The univariate structural time series model here is given by:
%     y_t = p_t + u_t + e_t
% p_{t+1} = p_t + b_t + c_t
% b_{t+1} = b_t + d_t
% u_t: trigonometric cycle
%
% The structural change in the trend occurs at tau = 2003.1, 
% but the change in slope already occurs at tau-1 = 2002.4:
%    p_{tau-1} =     p_{tau-2} + b_{tau-2} + c_{tau-2}
%    b_{tau-1} =         1 * w + b_{tau-2} + d_{tau-2}
%
%    p_tau =     p_{tau-1} + b_{tau-1} + c_{tau-1}
%    b_tau =                 b_{tau-1} + d_{tau-1}
%
% Given the structure of the system matrices W_t in the state space model, 
% it follows that:
%   W_{tau-2} = [0 1 0 0]', where tau is the time point of the intervention
%         W_t = [0 0 0 0]' for t ~= tau-2
% Construct super matrix W consisting of the matrices W_t

W = zeros(ly*4,1);
inW=in-2; 
W(((inW-1)*4+1):(((inW-1)*4)+1+3)) = [0 1 0 0]';

Y = [];                  % matrix of regression variables (stack of X_t
                         % matrices if X_t is time variant; X_t=X if X_t
                         % is time invariant.

npr=0;                   %number of forecasts
 
% Specify components, initial values and fixed parameters
comp.level=[1 0 0];
comp.slope=[1 0 0];
comp.irreg=[1  .1  NaN];
comp.cycle= [1 .1  NaN];
twopi=2*pi;
comp.cyclep=[0.9 twopi/40.; NaN NaN];
comp.cycleb=[twopi/60.  twopi/6.];
% comp.conout='cycle';    
comp.freq=freq;
comp.datei=datei;

% Put univariate structural time series model into state space form
[str,ferror] = suusm(comp,y,Y,npr);
% Add super matrix W (stack of W_t matrices for the state space model) to 
% the structre, str, of the state space model. By default, W=[]. 
str.W = W;

if ferror > 0
 return
end

% Estimate the state space model
[result,str] = usmestim(y,str);

%create estimated model
xvf=result.xvf;      % estimated parameters
xf=result.xf;        % fixed parameters
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result.tvr
%Note that the standard errors are divided by the concentrated parameter
%(sqrt(result.sigma2c))

%create estimated model
[X,Z,G,W,T,H,ins,ii,ferror] = pr2usm(xvf,xf,str);  

% %compute recursive residuals considering the regression coefficients fixed
% beta=result.h;        
% %beta is the vector containing the estimated regression coefficients. It 
% %will be considered fixed in the next function
% [Xtf,Ptf,recrs,recr,srecr,t1,A1,P1]=scakfff(y,X,Z,G,W,T,H,ins,ii,beta);
% %recrs contains the standardized residuals
% %t1 is the time at which the collapsed filter starts (no diffuse elements
% %in the Kalman filter)
% %A1 is x_{t1|t1-1}
% %P1 is MSE(A1)


% % Plot standardized recursive residuals (recrs)
% figure
% tsplot(recrs, datei, 'CWage(int_slope) St. Residuals');
% pause

%residual diagnostics
tname = 'CWage_STSM_sl_int_slope';
yor = y;
lam = 1;            %data are not transformed by taking logs
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
fname='results\CWage_STSM_sl_int_slope.txt';
fid = fopen(fname,'w');
%print estimation results
nreg=1;  
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
cyc = dispcomp(KKP,str,'cycle',datei,'GWage Smoothed Cycle');
pause
% saveas(gcf, 'results\CWage_STSM_sl_int_slope','fig');
% closefig

% Create vector of times
% bg_times = (bg_year+(bg_per/freq));
% end_times = (bg_year+(bg_per/freq)) + (ly/freq) - (1/freq);
% times = (bg_times:(1/freq):end_times)';

% Write the cycle along with the times vector into the Excel file
% xlswrite('results\CWage_STSM_int_slope.xls',[times cyc]);

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
