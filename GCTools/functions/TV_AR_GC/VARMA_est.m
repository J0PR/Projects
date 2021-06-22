function strv=VARMA_est(y,p,q)
if size(y,1)<size(y,2)%y must be a column vector
    y=y';
end
ndims=size(y,2);
morder= max(p,q);
kronIdxs=repmat(morder,1,ndims);
strv = estvarmaxkro(y,[],0,kronIdxs,0);

if p>q
    for i=1:p-q
        strv.theta(:,:,end-i+1)=0;
        strv.nparm=strv.nparm-ndims^2;
    end
elseif p<q
    for i=1:q-p
        strv.phi(:,:,end-i+1)=0;
        strv.nparm=strv.nparm-ndims^2;
    end
end

strv=mhanris(y,[],0,strv,0,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strv.lastOp='mhanris';
strv.lastPhi=strv.phis3;
strv.lastTheta=strv.thetas3;
strv.lastResid=strv.resid3;
strv.lastSigmar=strv.sigmar3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xvf,strv,ferror]=mconestim(y,[],strv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strv.lastOp='mconestim';
strv.lastPhi=strv.phiscon;
strv.lastTheta=strv.thetascon;
strv.lastResid=strv.residcon;
strv.lastSigmar=strv.sigmarcon;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=eye(ndims);
try
    [xvfx,strv,ferror]=mexactestimc(y,[],strv,Y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strv.lastOp='mexactestimc';
    strv.lastPhi=strv.phisexct;
    strv.lastTheta=strv.thetasexct;
    strv.lastResid=strv.e;
    strv.lastSigmar=strv.sigmarexct;
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
