function u=predError(X,alpha,nlags)
%get prediction error from data and MVAR coefficients

[nvar,nobs] = size(X);

alpha = alpha';
k=1;
for j=1:nlags
    for i=0:nvar-1
        beta(j+i*nlags,:)=alpha(k,:);
        k=k+1;
    end
end

% construct lag matrices
lags = -999*ones(nvar,nobs-nlags,nlags);
for jj=1:nvar
    for ii=1:nlags
        lags(jj,:,nlags-ii+1) = X(jj,ii:nobs-nlags+ii-1);
    end
end

%  regression (no constant term)
regressors = zeros(nobs-nlags,nvar*nlags);
for ii=1:nvar,
    s1 = (ii-1)*nlags+1;
    regressors(:,s1:s1+nlags-1) = squeeze(lags(ii,:,:));
end
for ii=1:nvar
    xvec = X(ii,:)';
    xdep = xvec(nlags+1:end);
    xpred(:,ii) = regressors*beta(:,ii);
    u(:,ii) = xdep-xpred(:,ii);
end