%Example of estimation of a transfer function model by exact maximum
%likelihood. The series is tf2, documented in SCA "Time Series for Research
%and Teaching".
%model is:
% (1-B)y_t = (3-2*B)(1-B)x_{t-1} + (1-.7*B)a_t
%
%input model is
% (1-B)x_t = alpha_t
%

%load data
yy=load('data\vf_tf2.dat'); y=yy(:,1); x=yy(:,2); 
seas=1;
[ny,s]=size(y);
[nx,mx]=size(x);

%difference series
yd=diferm(y,1); xd=diferm(x,1);

%estimate model using HR method (K.i. = 2)
kro=2; hr3=0; finv2=1;
strv = estvarmaxkro(yd,xd,seas,kro,hr3,finv2);


%fix unsignificant paramaters to zero and estimate again
strv.gamma(:,:,1)=0.; strv.phi(:,:,2:3)=zeros(1,2);
strv.theta(:,:,3)=0.;
strv.nparm=strv.nparm-4;
strv=mhanris(yd,xd,seas,strv,0,1); 


%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(yd,xd,strv);  
conp=strc.sigmarcon;

% %estimate model using the exact method
% Y=1;
% [xvfx,strx,ferror]=mexactestimc(yd,xd,strc,Y); 
% conp=strx.sigma2c;


%compute forecasts
npr=8; freq=1;
if (npr > 0)
 chb=1; Y=[];
 [ff,beta,e,f,strc,stx,recrs]=exactmedfvc(xvfc,yd,xd,strc,Y,chb);  
 %endogenous part
 A=stx.A; P=stx.P; Z=stx.Z; G=stx.G; T=stx.T; H=stx.H;
 hb=stx.hb; Mb=stx.Mb;
 Xp=Y;  
 Wp=[];
 cw=1.96;
 s=1;                           %number of series
 [pry,mypr,alpr,malpr]=ssmpred(npr,s,A,P,Xp,Z,G,Wp,T,H,hb,Mb);  
 spry=zeros(s,npr); 
 %exogenous part
 %inputs are stochastic
 hr3=1;  finv2=0;
 [strv,ferror] = estvarmaxpqrPQR(xd,[],freq,[0 0 0],[0 0 0],hr3,finv2);  
 sts.T=0; sts.Z=0; H=0; Sg=strv.sigmar2;
 [R,p]=chol(Sg); L=R'; sts.H=H; sts.G=L;
 [prx,mxpr,glpr,mglpr]=ssmpredexg(npr,xd,stx,sts);
 %forecasts and their mse
 pry=pry+prx; mypr=mypr*conp + mxpr;
 for i=1:npr
  spry(:,i)=sqrt(diag(mypr(:,:,i)));
 end
 opry=pry; ospry=spry;
 %plot forecasts 
 tname='tf2';
 out.pry=pry(1,:); out.spry=spry(1,:); 
 out.opry=opry(1,:); out.ospry=ospry(1,:); out.y=yd(:,1); 
 out.yor=yd(:,1); out.ny=length(yd(:,1)); out.npr=npr; out.cw=cw; 
 out.tname=tname;
 lam=1;                   %lam=0, logs are taken; =1, no logs are taken
 out.lam=lam; out.s=freq;
 pfctsusm(out);
end

