%Example of estimation of the flour price series using a
%VARMA model in echelon form for the differenced series. The model is 
%estimated by the conditional method.

%load data
y=load('data\flour-price.dat'); x=[];
%logs are taken
y=log(y); 
%series is differenced
yd=diferm(y,1);
seas=1;
[ny,s]=size(y);

lag=20; cw=1.96; ds=0; dr=0;
tname={'flour price 1','flour price 2','flour price 3'};
for i=1:3      
  c0=sacspacdif(yd(:,i),tname(i),dr,ds,seas,lag,cw);
  pause
end
closefig



%estimate model using HR method (K.i. = [1 0 0])
kro=[1 0 0];
strv = estvarmaxkro(yd,x,seas,kro,0,1); 

%check t-values
disp('t-values: ')
disp('tv-phi:')
disp(strv.phitv3)
disp('tv-theta:')
disp(strv.thetatv3)
disp('press any key to continue')
pause

%fix unsignificant paramaters to zero and estimate again
strv.phi(3,1,1)=0; strv.phi(1,1,2)=0; strv.theta(3,1,1)=0;
strv.nparm=strv.nparm-3;
strv=mhanris(yd,x,seas,strv,0,1);

%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(yd,x,strv);

%check t-values
disp('t-values: ')
disp('tv-phi:')
disp(strc.phitvcon)
disp('tv-theta:')
disp(strc.thetatvcon)
disp('press any key to continue')
pause

%compute autocovariance and autocorrelation matrices of residuals
lag=6; ic=1;
str=mautcov(strc.residcon,lag,ic);  
%see +-. signs
disp(str.sgn)
pause
disp(str.sgnt)
pause
