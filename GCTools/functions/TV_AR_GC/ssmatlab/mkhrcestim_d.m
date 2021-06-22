%Example of estimation of the simulated series used by Nsiri and Roy (1996)
%using a VARMA model in echelon form. The model is estimated by the 
%conditional method

clear
%load data
y=load('data\nsiri.dat'); x=[];
seas=1;
[ny,s]=size(y);

%estimate model using HR method (K.i. = [2 1])
strv = estvarmaxkro(y,x,seas,[2 1],0,1);  

%check t-values
disp('t-values: ')
disp('tv-phi:')
disp(strv.phitv3)
disp('tv-theta:')
disp(strv.thetatv3)
disp('press any key to continue')
pause


%fix unsignificant paramaters to zero and estimate again
strv.theta(2,:,2)=zeros(1,2);
strv.phi(1,:,3)=zeros(1,2); strv.nparm=strv.nparm-4;
strv=mhanris(y,x,seas,strv,0,1);

%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(y,x,strv);

%check t-values
disp('t-values: ')
disp('tv-phi:')
disp(strc.phitvcon)
disp('tv-theta:')
disp(strc.thetatvcon)
disp('press any key to continue')
pause

