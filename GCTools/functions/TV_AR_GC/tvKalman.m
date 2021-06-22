function ret=tvKalman(phi0,theta0,y)
% y(n*nPoints)
% n - Number of Vectors
% p - AR order
% q - MA order
% refs:
%   1) Analysis of earthquake ground motions using an improved Hilbert–Huang transform

[n n p]=size(phi0);

[n n q]=size(theta0);

[n nPoints]=size(y);

%fill state matrix: eq 12 Ref(1)
%|1 2 3 ... nPoints|
%|0 1 2 ... nPoints-1|
%x=zeros(n*p,nPoints+p-1);
x=zeros(n*p,nPoints);
for i=1:p
    x(1+(i-1)*n:i*n,i:end)=y(:,1:end-i+1);
end

v=zeros(n*q,nPoints+q-1);


I=eye(n);
phi=zeros([n n*p nPoints]);
theta=zeros([n n*q nPoints]);
A=zeros([n*p n*p nPoints]);
B=zeros([n*p n*q nPoints]);
C=zeros([n n*p nPoints]);
D=zeros([n n*q nPoints]);

phiArray=reshape(phi0,[n n*p]);
thetaArray=reshape(theta0,[n n*q]);

%initial states
X0=[phiArray thetaArray]';

%initial state posterior covariance
P0=10^4*eye(n*p+n*q);

%initial covariance matrix
R0=10^-4; %R0 can be set 10^-4 to 10^-2

%initial measurement vector
%H0=[y(1:p) zeros(1,q)];
H0=[x(:,p)' zeros(1,n*q)];

%covriance of process noise
Q=10^-4*eye(n*p+n*q);

ye=zeros(size(y));
e_pred_array=zeros(size(y));
errcov=zeros(1,nPoints);
PredErrCov=zeros([n n nPoints]);
Sigma_e=zeros([n n nPoints]);
for k=1:nPoints
    
    if k==1
        Ppost=P0;
        H=H0;
        Rpri=R0;
        Xpri=X0;
    else
        Xpri=Xpost;
        H=[x(:,k-1)' v(:,k-1)'];
        e_pred=y(:,k)-(H*Xpri)';
        %%%%%%%%%%teste%%%%%%%%%%%%%%
        ye(:,k) = H*Xpri;
        e_pred_array(:,k)=e_pred;
        if k>50%% for the covariance window <= 50
            ki=k-50;
        else
            ki=1;
        end
        PredErrCov(:,:,k)=cov(e_pred_array(:,ki:k)');
        errcov(k) = H*Ppost*H';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rpri=((k-1)*Rpost+e_pred'*e_pred)/k;
    end
    Ppri=Ppost+Q;
    K=Ppri*H'/(H*Ppri*H'+Rpri);
    Xpost = Xpri + K*(y(:,k)'-H*Xpri);

    Ppost = (eye(n*p+n*q)-K*H)*Ppri;
    e_resid=y(:,k)-(H*Xpost)';
    v=updateV(v,e_resid,k);
    %%%%%%%%%%TESTE2%%%%%%%%%%%%%
%     dowhile=1;
%     it=1;
%     while dowhile==1
%         if k==1
%             H=[x(:,k)' e_resid' zeros(1,n*(q-1))];
%         else
%             v(1:n,k-1)=e_resid;
%             H=[x(:,k-1)' v(:,k-1)'];%tenho de inserir o e_resid no v(1:n,k-1)
%         end
%         e_pred=y(:,k)-(H*Xpost)';
%         Rpri=((k-1)*Rpri+e_pred'*e_pred)/k;
%         Ppri=Ppost+Q;
%         K=Ppri*H'/(H*Ppri*H'+Rpri);
%         Xpost = Xpost + K*(y(:,k)'-H*Xpost);
%         Ppost = (eye(n*p+n*q)-K*H)*Ppri;
%         e_resid=y(:,k)-(H*Xpost)';
%         aux(:,it)=e_resid;
%         v=updateV(v,e_resid,k);
%         it=it+1;
%         if it>2000
%             dowhile=0;
%         end
%     end
    %%%%%%%%%%%FIMTESTE2%%%%%%%%%
    Rpost=Rpri;
    %%%%%%%%%%%FIMKalman%%%%%%%%%
    phi(:,:,k)=Xpost(1:n*p,1:n)';
    theta(:,:,k)=Xpost(n*p+1:end,1:n)';
    auxA1=eye(p-1);
    auxA2=[auxA1 zeros(p-1,1)];
    auxA=kron(auxA2,I);
    A(:,:,k)=[phi(:,:,k);auxA];
    
    auxB=zeros(p,1);
    auxB(1)=1;
    B(:,:,k)=kron(auxB,theta(:,:,k));
    
    auxC=zeros(1,p);
    auxC(1)=1;
    C(:,:,k)=kron(auxC,I);
    
    auxD=zeros(1,q);
    auxD(1)=1;
    D(:,:,k)=kron(auxD,I);
    
    
end
ret.A=A;
ret.B=B;
ret.C=C;
ret.D=D;
ret.phi=phi;
ret.theta=theta;
ret.ECov=PredErrCov;

function v=updateV(v,u,k)
n=size(u,1);
v(1:n,k)=u;
if k-1>0
    v(n+1:end,k)=v(1:end-n,k-1);
end


