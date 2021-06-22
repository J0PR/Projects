function [KKP,PT,recrs,recr,srecr,t1,A1,P1,KG]=scakfff(y,X,Z,G,W,T,H,ins,i,g)
%
%
%        This function applies the augmented Kalman filter and smoother 
%        to the series y corresponding to the model
%
%        y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%        alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%        collapse is applied to get rid of the diffuse component.
%
%        Input parameters:
%        y:     an (n x p) matrix of observations;
%        X    : an (n*p x nbeta) matrix containing the X_t matrices; 
%               a  (p x nbeta) if it is time invariant;
%               it can be []
%        Z    : an (n*p x nalpha) matrix containing the Z_t matrices;
%               a  (p x nalpha) matrix if it is time invariant
%        G    : an (n*p x nepsilon) matrix containing the G_t matrices;
%               a  (p x nepsilon) matrix if it is time invariant
%        W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%               an (nalpha x nbeta) matrix if it is time invariant;
%               it can be []
%        T    : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%               an (nalpha x nalpha) matrix if it time invariant
%        H    : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%               an (nalpha x nepsilon) if it is time invariant
%        ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%             state information, according to array i below
%        i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%             where
%             cc   = nalpha if c is not missing (0 if c missing)
%             cw0  = number of columns in W_0 (0 if W_0 missing)
%             ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%             cca1 = number of columns in A_1 (0 if A_1 missing)
%        g  : the beta vector, considered fixed
%
%        Output parameters:
%        KKP : an (n x nalpha) matrix containing the estimated x_{t|t}
%        PT  : an (n*nalpha x nalpha) matrix containing the
%              Mse of x_{t|t}
%       recrs: standardized recursive residuals
%        recr: recursive residuals
%       srecr: covariance matrices of recursive residuals
%          t1: initial time of the collapsed filter
%          A1: x_{t1|t1-1} initial state for the collapsed filter
%          P1: Mse(x_{t1|t1-1})
%          KG: stack of the Kalman gain vectors for the collapsed filter
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos, 
% Subdireccion Gral. de Analisis y P.E., 
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should 
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
%
%
% system dimensions
%
[n, p]=size(y);
[mx, nx]=size(X);
[mz, nalpha]=size(Z);
[mg, neps]=size(G);
[mw, nw]=size(W);
[mt, nt]=size(T);
[mh, nh]=size(H);
[mi, ni]=size(ins);
[mc,nc]=size(i);
nbeta=max(nx,nw);
cc=i(1);
cw0=i(2);
ca1=i(3);
cca1=i(4);
%
% check for inconsistencies
%
if mi ~= nalpha
 disp('the number of rows in ins is incorrect')
 return
end
if mx > p || mz > p || mg > p || mw > nalpha
   ...| mt > nalpha | mh > nalpha
%
% system matrices are time varying
%
   if mx > p && mx ~= n*p
      disp('the number of rows in X is incorrect')
      return
   end
   if mz > p && mz ~= n*p
      disp('the number of rows in Z is incorrect')
      return
   end
   if mg > p && mg ~= n*p
      disp('the number of rows in G is incorrect')
      return
   end
   if mw > nalpha && mw ~= n*nalpha
      disp('the number of rows in W is incorrect')
      return
   end
   if mt > nalpha && mt ~= n*nalpha
      disp('the number of rows in T is incorrect')
      return
   end
   if mh > nalpha && mh ~= n*nalpha
      disp('the number of rows in H is incorrect')
      return
   end
end
if nw ~= nx && nw ~= 0 && nx ~= 0
  disp('the number of columns in W and X should agree')
  return
end
if nt ~= nalpha
  disp('the number of columns in T and Z should agree')
  return
end
if nh ~= neps
  disp('the number of columns in H and G should agree')
  return
end
if mi ~= nalpha
  disp('the number of rows in ins and the number of columns in')
  ... disp(' Z should agree')
  return
end
if nc ~= 4 || mc ~= 1
  disp('c should be a 1 x 4 matrix')
  return
end
if nbeta == 0 && cw0 ~= 0
  disp('the number of rows in W_0 should be zero')
  return
end
if ni ~= cc+cw0+ca1+cca1
  disp('the number of columns in ins should be ')
  disp(cc+cw0+ca1+cca1)
  return
end
%
% initial covariance matrix
%
if cc == 0
   P=zeros(nalpha,nalpha);
else
   P=ins(:,1:cc);
end
%
% initial state
%
if cw0 == 0
   if nbeta == 0
    W0=[];
   else
    W0=zeros(nalpha,nbeta);
   end
else
   W0=ins(:,cc+1:cc+cw0);
end
%W0
if ca1 == 0
   a1=zeros(nalpha,1);
else
   a1=ins(:,cc+cw0+1:cc+cw0+ca1);
end
if cca1 == 0
   ndelta=0;
   aa1=[];
else
   ndelta=cca1;
   aa1=ins(:,cc+cw0+ca1+1:cc+cw0+ca1+cca1);
end
A=[-aa1 -full(W0) a1];
%
% augmented Kalman filter
%
%A
%P
nd=ndelta;
ndb=ndelta+nbeta;
ndb1=ndb+1;
nb1=nbeta+1;
KKP=zeros(n,nalpha);
PT=zeros(n*nalpha,nalpha);
SQT=zeros(ndb+p,ndb1);
PT(1:nalpha,:)=P;
%tol=1.d-10;
iti=0;
collps=0;
%the following line added on 5-6-2010
if nd == 0, collps=1; end
%end of addition
%the following line modified 1.2012
qt=0; st2=0; ndbc=nd; nmiss=0;
%the following line added 1-2012
recrs=[]; recr=[]; srecr=[]; KG=[];
%the following three lines added 1-2012
if collps == 1
 A1=A*[-g' 1]'; P1=P;
end
%
   for i=1:n
%   i
    ip=(i-1)*p+1:i*p;
    ia=(i-1)*nalpha+1:i*nalpha;
    if mx > p 
       XX=X(ip,:);
    else
       XX=X;
       if nbeta > 0 && isempty(X) == 1
        XX=zeros(p,nbeta);
       end
    end
    if mz > p 
       ZZ=Z(ip,:);
    else
       ZZ=Z;
    end
    if mg > p 
       GG=G(ip,:);
    else
       GG=G;
    end
    if mw > nalpha 
       WW=W(ia,:);
    else
       WW=W;
       if nbeta > 0 && isempty(W) == 1
        WW=zeros(nalpha,nbeta);
       end
    end
    if mt > nalpha
       TT=T(ia,:);
    else
       TT=T;
    end
    if mh > nalpha
       HH=H(ia,:);
    else
       HH=H;
    end
    YY=y(i,:)';
%
% check for missing values
%
    [J,miss]=findmv(YY);
    nmiss=nmiss+miss;
    if miss > 0 && miss < p
     J=sparse(J);
     if ~isempty(XX)
      XX=J*XX;
     end
     YY=J*YY;
     if ~isempty(GG)
      GG=J*GG;
     end
     if ~isempty(ZZ)
      ZZ=J*ZZ;
     end
     V=[zeros(p-miss,ndelta) full(XX) YY] - ZZ*A;
     F=ZZ*P*ZZ' + GG*GG'; [CF,cp]=chol(F);
     if cp == 0
      Frt=CF';
%       Lm1=inv(Frt);
%       DD=Lm1*V;
%       K=(TT*P*ZZ'+HH*GG')*Lm1'*Lm1;
      DD=Frt\V;
      K=(TT*P*ZZ'+HH*GG')/F;
     elseif cp == 1
%modified 19th Nov. 2004      
%       DD=zeros(p,ndb1);
      if ndelta == 0, DD=[]; else DD=zeros(p,ndb1); end
%end of modification      
      K=zeros(nalpha,p-miss);Lm1=zeros(p); %last statement added on 18th Dec. 2005
     else
      error('singular matrix different from zero in scakfff')
     end
    elseif miss == p
%modified 19th Nov. 2004      
%       DD=zeros(p,ndb1);
      if ndelta == 0, DD=[]; else DD=zeros(p,ndb1); end
%end of modification      
     K=zeros(nalpha,p);Lm1=zeros(p); %last statement added on 18th Dec. 2005
    else
     V=[zeros(p,ndelta) full(XX) YY] - ZZ*A;
     F=ZZ*P*ZZ' + GG*GG'; [CF,cp]=chol(F);
     if cp == 0
      Frt=CF';
%       Lm1=inv(Frt);
%       DD=Lm1*V;
%       K=(TT*P*ZZ'+HH*GG')*Lm1'*Lm1;
      DD=Frt\V;
      K=(TT*P*ZZ'+HH*GG')/F;
     elseif cp == 1
%modified 19th Nov. 2004      
%       DD=zeros(p,ndb1);
      if ndelta == 0, DD=[]; else DD=zeros(p,ndb1); end
%end of modification      
      K=zeros(nalpha,p);Lm1=zeros(p); %last statement added on 18th Dec. 2005
     else
      error('singular matrix different from zero in scakfff')    
     end
    end
%
% filtering
%
    ia=(i-1)*nalpha+1:i*nalpha;
    if (nd == 0) || (collps == 1 && i > iti)
      if  miss < p  
%       U=SQT(1:nbeta,1:nbeta); 
%       invu=inv(U); 
%       deno=i*p-nmiss-ndbc; % if deno > 10*ndbc, deno=deno-ndbc; end
      deno=i*p-nmiss-ndbc;
      et=V*[-g' 1]'; 
      set=F;
%the following three lines added 1-2012      
      ets=Frt\et;  recrs=[recrs; ets']; 
      recr=[recr; et']; srecr=[srecr; set];
      KG=[KG; K];
%end of addition
      qt=qt+et'*(set\et); 
%       qt=qt+et'*et; 
      st2=qt/deno; 
      RNN=(Frt\ZZ)';
      R=RNN*DD;  NN=RNN*RNN';
      KKP(i,:)=((A+P*R)*[-g' 1]')';
%       DB=A(:,1:ndb)+P*R(:,1:ndb);
      PT(ia,:)=(P-P*NN*P)*st2;
     elseif  miss == p  
      KKP(i,:)=(A*[-g' 1]')';
      PT(ia,:)=P*st2;
     else
      KKP(i,:)=NaN;
      PT(ia,:)=NaN;
     end
    else
     KKP(i,:)=NaN;
     PT(ia,:)=NaN;
    end
    if miss == p
     A=[zeros(nalpha,ndelta) -full(WW) zeros(nalpha,1)]+TT*A;
     P=TT*P*TT'+HH*HH';
    else
     A=[zeros(nalpha,ndelta) -full(WW) zeros(nalpha,1)]+TT*A+K*V;
     P=TT*P*(TT-K*ZZ)'+HH*(HH-K*GG)';
    end
%
% SQT updating
%
% the following line changed on 1-2012
    if (ndb > 0) && (collps == 0)
     [Q,SQT]=qr([DD;SQT(1:ndb,:)]); 
%      if nbeta > 0 & best == 0 & collps == 1
%         if min(abs(diag(SQT(ndelta+1:ndb,ndelta+1:ndb)))) > tol 
%          best = 1; initf = i+1;
%         end
%      end
    else
     SQT=DD;
    end
%
% single collapse if possible
%
    if collps == 0
%     if min(abs(diag(SQT(1:ndelta,1:ndelta)))) > tol
      if rank(SQT(1:ndelta,1:ndelta)) == ndelta
%     i
%     SQT
%     Q
      Rm1=pinv(SQT(1:nd,1:nd));
      hd=Rm1*SQT(1:nd,nd+1:ndb1);
      M=Rm1*Rm1';
      AA=A(:,1:nd);
      P=P+AA*M*AA';
      A=A(:,nd+1:ndb1) - AA*hd;
      %the following line added 1-2012
      A1=A*[-g' 1]'; P1=P;
%
% redefinition of SQT
%
%     disp('SQT after collapsing')
      [nsqt,msqt]=size(SQT);
      SQT=SQT(nd+1:nsqt,nd+1:ndb1);
      ndelta=0;
      ndb1=nb1;
      ndb=nbeta;
      iti=i;
      collps=1;
     end
    end
%   pause
   end
   t1=iti+1;
%SQT
%KKP
%PT
