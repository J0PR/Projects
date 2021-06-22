function [ret] = mvar_0lag_in(X,nlags,opt,STATFLAG)
X=X';
if nargin<4,
    STATFLAG = 1;
end
if nargin<3,
    opt='LS';
end
if strcmp(opt,'CP')
    nlagsP1=nlags+1;
else
    nlagsP1=nlags;
end
% figure regression parameters
[nvar,nobs] = size(X);
if(nvar>nobs) error('nvar>nobs, check input matrix'); end

% remove sample means if present (no constant terms in this regression)
m = mean(X');
if(abs(sum(m)) > 0.0001)
    mall = repmat(m',1,nobs);
    X = X-mall;
end

% construct lag matrices
lags = -999*ones(nvar,nobs-nlagsP1,nlagsP1);
for jj=1:nvar
    for ii=1:nlagsP1
        lags(jj,:,nlagsP1-ii+1) = X(jj,ii:nobs-nlagsP1+ii-1);
    end
end

%  regression (no constant term)
regressors = zeros(nobs-nlagsP1,nvar*nlagsP1);
for ii=1:nvar,
    s1 = (ii-1)*nlagsP1+1;
    regressors(:,s1:s1+nlagsP1-1) = squeeze(lags(ii,:,:));
end
if strcmp(opt, 'CP')
    new_beta=zeros(nvar*nlagsP1,nvar);
    for ii=1:nvar
        xvec = X(ii,:)';
        xdep = xvec(nlags+1:end-1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_regressors=regressors;
        loc=(ii-1)*nlagsP1+1;
        temp_regressors(:,loc)=[];
        beta(:,ii) = temp_regressors\xdep;
        
        new_beta(1:loc-1,ii)=beta(1:loc-1,ii);
        new_beta(loc,ii)=0;
        new_beta(loc+1:end,ii)=beta(loc:end,ii);
        xpred(:,ii) = temp_regressors*beta(:,ii);  % old beta to have the same number of vars as tem_regressors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     beta(:,ii) = regressors\xdep;
        %     xpred(:,ii) = regressors*beta(:,ii);  % keep hold of predicted values
        u(:,ii) = xdep-xpred(:,ii);
        RSS(ii) = sum(u(:,ii).^2);        
    end
    beta = new_beta;
    k=1;
    for j=1:nlagsP1
        for i=0:nvar-1
            alpha(k,:)=beta(j+i*nlagsP1,:);
            k=k+1;
        end
    end
    alpha = alpha';
%     alpha(:,1:nvar)=[];
%     beta(1:nlagsP1:nvar*nlagsP1)=[];
else
    for ii=1:nvar
        xvec = X(ii,:)';
        xdep = xvec(nlags+1:end);
        beta(:,ii) = regressors\xdep;
        xpred(:,ii) = regressors*beta(:,ii);  % keep hold of predicted values
        u(:,ii) = xdep-xpred(:,ii);
        RSS(ii) = sum(u(:,ii).^2);
    end
    k=1;
    for j=1:nlags
        for i=0:nvar-1
            alpha(k,:)=beta(j+i*nlags,:);
            k=k+1;
        end
    end
    alpha = alpha';
end
% check whiteness & consistency if required
if STATFLAG,
    for ii=1:nvar,
        waut(ii) = cca_whiteness(X,u(:,ii));
    end
    cons = cca_consistency(X,xpred);
else
    waut = -1;
    cons = -1;
end

%   organize output structure
ret.beta = beta;
ret.alpha=alpha;
ret.u = u;
ret.rss = RSS;
ret.cons = cons;
ret.waut = waut;
ret.Z = cov(u);

%    Written by João Rodrigues
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


