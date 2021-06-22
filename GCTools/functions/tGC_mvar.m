function [ret] = tGC_mvar(X,nlags,extended_AR)
% -----------------------------------------------------------------------
%   FUNCTION: cca_granger_regress.m
%   PURPOSE:  perform multivariate regression with granger causality statistics
%
%   INPUT:  X           -   nvar (rows) by nobs (cols) observation matrix
%           nlags       -   number of lags to include in model
%           STATFLAG    -   if 1 do F-stats
%
%   OUTPUT: ret.covu    -   covariance of residuals for unrestricted model
%           ret.covr    -   covariance of residuals for restricted models
%           ret.prb     -   Granger causality probabilities (column causes
%                           row. NaN for non-calculated entries)
%           ret.fs      -   F-statistics for above
%           ret.gc      -   log ratio causality magnitude
%           ret.doi     -   difference-of-influence (based on ret.gc)
%           ret.rss     -   residual sum-square
%           ret.rss_adj -   adjusted residual sum-square
%           ret.waut    -   autocorrelations in residuals (by Durbin
%           Watson)
%           ret.cons    -   consistency check (see cca_consistency.m)

%   Written by Anil K Seth Sep 13 2004
%   Updated AKS December 2005
%   Updated AKS November 2006
%   Updated AKS December 2007 to do ratio based causality
%   Updated AKS May 2008, fix nlags = 1 bug.
%   Updated AKS Apr 2009, difference-of-influence and optional stats
%   Updated AKS Aug 2009, specify regressor matrix size in advance
%   Updated AKS Aug 2009, implement whiteness + consistency checks
%   Ref: Seth, A.K. (2005) Network: Comp. Neural. Sys. 16(1):35-55
% COPYRIGHT NOTICE AT BOTTOM
% -----------------------------------------------------------------------

% SEE COPYRIGHT/LICENSE NOTICE AT BOTTOM

% figure regression parameters
nobs = size(X,2);
nvar = size(X,1);
if(nvar>nobs) error('error in cca_granger_regress: nvar>nobs, check input matrix'); end
STATFLAG = 1;


[ARF,RCF,PE,DC,u] = mvar(X', nlags);
if extended_AR
    M = size(ARF,1);
    covM  = PE(:,(1-M:0)+end);
    [L,D]=ldl(covM);
    u=(L\u')';
end
for ii=1:nvar
    C(ii) = covariance(u(:,ii),u(:,ii),nobs-nlags);
    RSS1(ii) = sum(u(:,ii).^2);
end
covu = cov(u);

%   A rectangular matrix A is rank deficient if it does not have linearly independent columns.
%   If A is rank deficient, the least squares solution to AX = B is not unique.
%   The backslash operator, A\B, issues a warning if A is rank deficient and
%   produces a least squares solution that has at most rank(A) nonzeros.

%   restricted regressions (no constant terms)
for ii=1:nvar
    xvec = X(ii,:)';
    xdep = xvec(nlags+1:end);          % dependent variable
    caus_inx = setdiff(1:nvar,ii);     % possible causal influences on xvec
    u_r = zeros(nobs-nlags,nvar,'single');
    for jj=1:length(caus_inx)
        eq_inx = setdiff(1:nvar,caus_inx(jj));  % vars to include in restricted regression (jj on ii)
        [ARF,RCF,PE,DC,temp_temp_r] = mvar(X(eq_inx,:)', nlags);
        if extended_AR
            M = size(ARF,1);
            covM  = PE(:,(1-M:0)+end);
            [L,D]=ldl(covM);
            temp_temp_r=(L\temp_temp_r')';
        end
        temp_r=temp_temp_r(:,find(eq_inx==ii));
        RSS0(ii,caus_inx(jj)) = sum(temp_r.^2);
        S(ii,caus_inx(jj)) = covariance(temp_r,temp_r,nobs-nlags); % dec 08
        u_r(:,caus_inx(jj)) = temp_r;        
    end
    covr{ii} = cov(u_r);
end

% calc Granger values
gc = ones(nvar).*NaN;
igc = ones(nvar).*NaN;
doi = ones(nvar).*NaN;
%   do Granger f-tests if required
if STATFLAG == 1,
    prb = ones(nvar).*NaN;
    ftest = zeros(nvar);
    n2 = (nobs-nlags)-(nvar*nlags);
    for ii=1:nvar-1
        for jj=ii+1:nvar
            ftest(ii,jj) = ((RSS0(ii,jj)-RSS1(ii))/nlags)/(RSS1(ii)/n2);    % causality jj->ii
            prb(ii,jj) = 1 - cca_cdff(ftest(ii,jj),nlags,n2);
            ftest(jj,ii) = ((RSS0(jj,ii)-RSS1(jj))/nlags)/(RSS1(jj)/n2);    % causality ii->jj
            prb(jj,ii) = 1 - cca_cdff(ftest(jj,ii),nlags,n2);
            gc(ii,jj) = log(S(ii,jj)/C(ii));
            gc(jj,ii) = log(S(jj,ii)/C(jj));
            doi(ii,jj) = gc(ii,jj) - gc(jj,ii);
            doi(jj,ii) = gc(jj,ii) - gc(ii,jj);
            igc(ii,jj) = log(S(ii,jj)*S(jj,ii)/det(covu));
            igc(jj,ii)=igc(ii,jj);
        end
    end
else
    ftest = -1;
    prb = -1;
    for ii=1:nvar-1,
        for jj=ii+1:nvar,
            gc(ii,jj) = log(S(ii,jj)/C(ii));
            gc(jj,ii) = log(S(jj,ii)/C(jj));
            igc(ii,jj) = log(S(ii,jj)*S(jj,ii)/det(C));
            igc(jj,ii)=igc(ii,jj);
            doi(ii,jj) = gc(ii,jj) - gc(jj,ii);
            doi(jj,ii) = gc(jj,ii) - gc(ii,jj);
        end
    end
end

%   do r-squared and check whiteness, consistency

%   organize output structure
ret.gc = gc;
ret.igc = igc;
ret.fs = ftest;
ret.prb = prb;
ret.covu = covu;
ret.covr = covr;
ret.doi = doi;
ret.type = 'td_normal';

% This file is part of GCCAtoolbox.  It is Copyright (C) Anil Seth, 2004-09
% 
% GCCAtoolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GCCAtoolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GCCAtoolbox.  If not, see <http://www.gnu.org/licenses/>.

