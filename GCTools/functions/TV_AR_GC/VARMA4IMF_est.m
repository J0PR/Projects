function strv=VARMA4IMF_est(y,p,q,nIMFs)
ndims=size(y,2);
nVars= ndims/nIMFs;
morder= max(p,q);
kronIdxs=repmat(morder,1,ndims);
strv = estvarmaxkro(y,[],0,kronIdxs,0);

if p>q
    for i=1:p-q
        strv.theta(:,:,end)=0;
        strv.nparm=strv.nparm-ndims^2;
    end
elseif p<q
    for i=1:q-p
        strv.phi(:,:,end)=0;
        strv.nparm=strv.nparm-ndims^2;
    end
end
%%%%%%%% for IMFs modification: only self predict %%%%%%%%
eyeNaN=zeros(nIMFs);
eyeNaN(find(eye(nIMFs)))=NaN;
MNaN=kron(ones(nVars),eyeNaN);
%os phi e theta tem dim 3 = p+1 e q+1 respectivamente pois têm sempre uma
%matriz identidade no inicio (lag 0).
strv.phi(:,:,2:1+p)=repmat(MNaN,[1 1 p]);
strv.theta(:,:,2:1+q)=repmat(MNaN,[1 1 q]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strv=mhanris(y,[],0,strv,0,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strv.lastOp='mhanris';
strv.lastPhi=strv.phis3;
strv.lastTheta=strv.thetas3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xvf,strv,ferror]=mconestim(y,[],strv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strv.lastOp='mconestim';
strv.lastPhi=strv.phiscon;
strv.lastTheta=strv.thetascon;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=eye(ndims);
try
    [xvfx,strv,ferror]=mexactestimc(y,[],strv,Y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strv.lastOp='mexactestimc';
    strv.lastPhi=strv.phisexct;
    strv.lastTheta=strv.thetasexct;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catch err
    disp('Error in: [xvfx,strv,ferror]=mexactestimc(y,[],strv,Y);')
    disp(err.message)
end
%[xvfx,strv,ferror]=mexactestimc(y,[],strv,Y);

if p>q
    for i=1:p-q
        strv.lastTheta(:,:,end)=[];
    end
elseif p<q
    for i=1:q-p
        strv.lastPhi(:,:,end)=[];
    end
end
