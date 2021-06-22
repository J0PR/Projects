function [ARF,RCF,PE,DC,F] = mvar(Y, Pmax)
% MVAR estimates Multi-Variate AutoRegressive model parameters
% Partial Correlation Estimation: Vieira-Morf Method with unbiased covariance estimation
% all estimators can handle data with missing values encoded as NaNs.
%
% 	[AR,RC,PE] = mvar(Y, p);
% 	[AR,RC,PE] = mvar(Y, p, Mode);
%
% INPUT:
%  Y	 Multivariate data series
%  p     Model order
%
% OUTPUT:
%  AR    multivariate autoregressive model parameter
%  RC    reflection coefficients (= -PARCOR coefficients)
%  PE    remaining error variance
%  F     forward error
%
% All input and output parameters are organized in columns, one column
% corresponds to the parameters of one channel.
%
%
%
% A multivariate inverse filter can be realized with
%   [AR,RC,PE] = mvar(Y,P);
%   e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR,1)),Y);
%
% see also: MVFILTER, MVFREQZ, COVM, SUMSKIPNAN, ARFIT2

%	$Id: mvar.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2006 by Alois Schloegl <a.schloegl@ieee.org>
%       This is part of the TSA-toolbox. See also
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/
%
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



[N,M] = size(Y);

if nargin<2,
    Pmax=max([N,M])-1;
end;

if iscell(Y)
    Pmax = min(max(N ,M ),Pmax);
    C    = Y;
end;

[C(:,1:M),n] = covm(Y,'M');
PE(:,1:M)  = C(:,1:M)./n;


F = Y;
B = Y;
PEF = PE(:,1:M);
PEB = PE(:,1:M);
for K = 1:Pmax,
    [D,n]	= covm(F(K+1:N,:),B(1:N-K,:),'M');
    D = D./n;
    
    ARF(:,K*M+(1-M:0)) = D / PEB;
    ARB(:,K*M+(1-M:0)) = D'/ PEF;
    
    tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
    B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
    F(K+1:N,:) = tmp;
    
    for L = 1:K-1,
        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
        ARF(:,L*M+(1-M:0))   = tmp;
    end;
    
    RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
    RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
    
    [PEF,n] = covm(F(K+1:N,:),F(K+1:N,:),'M');
    PEF = PEF./n;
    
    [PEB,n] = covm(B(1:N-K,:),B(1:N-K,:),'M');
    PEB = PEB./n;
    
    PE(:,K*M+(1:M)) = PEF;
end;
F=F(Pmax+1:end,:);

if any(ARF(:)==inf),
    p = 3;
    FLAG_MATRIX_DIVISION_ERROR = ~all(all(isnan(repmat(0,p)/repmat(0,p)))) | ~all(all(isnan(repmat(nan,p)/repmat(nan,p))));
    
    if FLAG_MATRIX_DIVISION_ERROR,
        %fprintf(2,'### Warning MVAR: Bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.\n');
        warning('MVAR: bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.');
        
        %%%%% Workaround
        ARF(ARF==inf)=NaN;
        RCF(RCF==inf)=NaN;
    end;
end;

%MAR   = zeros(M,M*Pmax);
DC     = zeros(M);
for K  = 1:Pmax,
    %       VAR{K+1} = -ARF(:,K*M+(1-M:0))';
    DC  = DC + ARF(:,K*M+(1-M:0))'.^2; %DC meausure [3]
end;

function [CC,NN] = covm(X,Y,Mode)
% COVM generates covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
%
% COVM(X,Mode);
%      calculates the (auto-)correlation matrix of X
% COVM(X,Y,Mode);
%      calculates the crosscorrelation between X and Y
%
% Mode = 'M' minimum or standard mode [default]
% 	C = X'*X; or X'*Y correlation matrix
%
% Mode = 'E' extended mode
% 	C = [1 X]'*[1 X]; % l is a matching column of 1's
% 	C is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
% 	the mean (or sum) is stored on the 1st row and column of C
%
% Mode = 'D' or 'D0' detrended mode
%	the mean of X (and Y) is removed. If combined with extended mode (Mode='DE'), 
% 	the mean (or sum) is stored in the 1st row and column of C. 
% 	The default scaling is factor (N-1). 
% Mode = 'D1' is the same as 'D' but uses N for scaling. 
%
% C = covm(...); 
% 	C is the scaled by N in Mode M and by (N-1) in mode D.
% [C,N] = covm(...);
%	C is not scaled, provides the scaling factor N  
%	C./N gives the scaled version. 


if nargin<3,
        if nargin==2,
		if isnumeric(Y),
                        Mode='M';
                else
		        Mode=Y;
                        Y=[];
                end;
        elseif nargin==1,
                Mode = 'M';
                Y = [];
        elseif nargin==0,
                error('Missing argument(s)');
        end;
end;        

Mode = upper(Mode);

[r1,c1]=size(X);
if ~isempty(Y)
        [r2,c2]=size(Y);
        if r1~=r2,
                error('X and Y must have the same number of observations (rows).');
                return;
        end;
else
        [r2,c2]=size(X);
end;

if (c1>r1) | (c2>r2),
        warning('Covariance is ill-defined, because of too less observations (rows)');
end;

if ~isempty(Y),
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	NN = real(~isnan(X)')*real(~isnan(Y));
	        X(isnan(X)) = 0; % skip NaN's
	        Y(isnan(Y)) = 0; % skip NaN's
        	CC = X'*Y;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S1,N1] = sumskipnan(X,1);
                [S2,N2] = sumskipnan(Y,1);
                
                NN = real(~isnan(X)')*real(~isnan(Y));
        
	        if any(Mode=='D'), % detrending mode
        		X  = X - ones(r1,1)*(S1./N1);
                        Y  = Y - ones(r1,1)*(S2./N2);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else   %  'D0'       
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(isnan(X)) = 0; % skip NaN's
        	Y(isnan(Y)) = 0; % skip NaN's
                CC = X'*Y;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N2; N1', NN];
                        CC = [r1, S2; S1', CC];
                end;
	end;
        
else        
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	tmp = real(~isnan(X));
                NN  = tmp'*tmp; 
                X(isnan(X)) = 0; % skip NaN's
	        CC = X'*X;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S,N] = sumskipnan(X,1);
        	tmp = real(~isnan(X));
                NN  = tmp'*tmp; 
                if any(Mode=='D'), % detrending mode
	                X  = X - ones(r1,1)*(S./N);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else  %  'D0'      
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(isnan(X)) = 0; % skip NaN's
                CC = X'*X;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N; N', NN];
                        CC = [r1, S; S', CC];
                end;
	end
end;

if nargout<2
        CC = CC./NN; % unbiased
end;
return; 