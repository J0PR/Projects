%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.6 of Reinsel (1997), pp. 221-224
%
% Series are: 1) US housing starts and 2) US housing sold.
% The period is: January 1965 through December 1974.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

yt=load('c:\ssm_matlab\data\housing.dat'); 
y=yt(:,3:4);

lag=36; cw=1.96; freq=12; dr=0; 
tname={'US housing starts','US housing sold'};
for i=1:2      
 for ds=0:1
  c0=sacspacdif(y(:,i),tname(i),dr,ds,freq,lag,cw);
  pause
 end
end
closefig


yd=diferm(y,freq);

%estimate a VARMAX(1,0,0)(0,1,0) model by the Hannan-Rissanen method. 
x=[]; hr3=0; finv2=1;  
[strv,ferror] = estvarmaxpqrPQR(yd,x,freq,[1 0 0],[0 1 0],hr3,finv2); 

disp('estimated phi parameters:')
disp(strv.phis3(:,:,1:2))
pause
disp('estimated theta parameters:')
disp(strv.thetas3(:,:,1:13))
pause


%estimate using exact ML
%setup model
Phi=eye(2); Th(:,:,1)=eye(2); Th(:,:,2)=strv.thetas3(:,:,13);
phi=strv.phis3(:,:,1:2); th=eye(2); Sigma=strv.sigmar3;  

[str,ferror] = suvarmapqPQ(phi,th,Phi,Th,Sigma,freq); 

Y=[]; 

%estimate model using the exact method
result=varmapqPQestim(yd,str,Y);  

%estimated and fixed parameters
xvf=result.xvf; xf=result.xf;  
%t-values of varma estimated parameters are in result.tv

%create estimated model
[phif,thf,Phif,Thf,Lf,ferror] = pr2varmapqPQ(xvf,xf,str);
Sigmar=result.Sigmar;
%t-values
tvf=result.tv; 
[phitvf,thtvf,Phitvf,Thtvf,Ltvf,ferror] = pr2varmapqPQ(tvf,xf,str);  

disp('estimated phi:')
disp(phif)
pause
disp('estimated Th:')
disp(Thf)
pause

disp('t-values of phi:')
disp(phitvf)
pause
disp('t-values of Th:')
disp(Thtvf)
pause

disp('Estimated covariance matrix of residuals:')
disp(result.Sigmar)


