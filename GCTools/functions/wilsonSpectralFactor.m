function [H, Z, S, psi] = wilsonSpectralFactor(S) %#codegen
% Usage  : [H, Z, S, psi] = wilsonSpectralFactor(S,freq);
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : psi (left spectral factor)
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Addapted from M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.
tol=1e-18;

% number of channels
m   = size(S,1);
N   = size(S,3)-1;
N2  = 2*N;

% preallocate memory for efficiency
Sarr   = complex(zeros(m,m,N2));
gam    = zeros(m,m,N2);
gamtmp = zeros(m,m,N2);
psi    = complex(zeros(m,m,N2));
psi_alt    = complex(zeros(m,m,N2));
nanfree_psi= complex(zeros(m,m,N2));
I      = eye(m); % Defining m x m identity matrix
Niterations=200;
%Step 1: Forming 2-sided spectral densities for ifft routine in matlab

%for f = freq
for f_ind = 1:size(S,3)
  Sarr(:,:,f_ind) = S(:,:,f_ind);
  if(f_ind>1)
    Sarr(:,:,2*N+2-f_ind) = S(:,:,f_ind).';
  end
end

%Step 2: Computing covariance matrices
for k1 = 1:m
  for k2 = 1:m
    gam(k1,k2,:) = real(ifft(squeeze(Sarr(k1,k2,:))));
  end
end

%Step 3: Initializing for iterations 
% gam0 = gam(:,:,1);
% h = my_chol(gam0);
% [h, dum] = chol(gam0);
% if dum
%   h = rand(m,m); h = triu(h); %arbitrary initial condition
% end
h = rand(m,m); h = triu(h);
for ind = 1:N2
  psi(:,:,ind) = h; 
end
psierrfVect=zeros(1,Niterations);
toterrfVect=zeros(1,Niterations);
%Step 4: Iterating to get spectral factors
err_alt=0;
alt=0;
iter=0;
err_best=Inf;
psi_best=[];
while iter <Niterations
    iter=iter+1;
    g = complex(zeros(size(psi)));
  for ind = 1:N2
    %invpsi     = pinv(psi(:,:,ind));% + I*eps(psi(:,:,ind))); 
    %g(:,:,ind) = invpsi*Sarr(:,:,ind)*invpsi'+I;%Eq 3.1
    g(:,:,ind) = psi(:,:,ind)\Sarr(:,:,ind)/psi(:,:,ind)'+I;%Eq 3.1,
  end
  gp = PlusOperator(g,m,N+1); %gp constitutes positive and half of zero lags 

  psi_old = psi;
  if sum(isnan(psi(:))) == 0
      nanfree_psi = psi;
  end
  psierr=zeros(1,N2);
  toterr=zeros(1,N2);
  for k = 1:N2
    psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
    psierr(k)  = norm(psi(:,:,k)-psi_old(:,:,k),1);
    toterr(k)  = norm(psi(:,:,k)*psi(:,:,k)'-Sarr(:,:,k),1);
  end
  psierrf = mean(psierr);
  toterrf = mean(toterr);
  psierrfVect(iter)=psierrf;
  toterrfVect(iter)=toterrf;
  
  if iter ==1
      psi_best=psi;
  end
  
  if toterrf<=err_best
      psi_best=psi;
      err_best=toterrf;
  end
  
  if iter>1
      if toterrfVect(iter)>10*toterrfVect(iter-1)
          if iter<25
              if alt==0
                  alt=1;
                  err_alt=toterrfVect(iter-1);
              else
                  if err_alt>=toterrfVect(iter-1)
                      psi=psi_old;
                      break;
                  else
                      psi=psi_alt;
                      break;
                  end
              end
              h = rand(m,m); h = triu(h);
              psi_alt=psi_old;
              for ind = 1:N2
                  psi(:,:,ind) = h;
              end
              iter=0;
              continue;
          else
              if alt~=0
                  if err_alt>=toterrfVect(iter-1)
                      psi=psi_old;
                      break;
                  else
                      psi=psi_alt;
                      break;
                  end
              end
              psi=psi_old;
              break;
          end
      end
  end
  if iter>25
      if toterrfVect(iter)>0.001
          if abs(toterrfVect(iter)-toterrfVect(iter-3))<0.00001
              break;
          end
      end
  end
  if(psierrf<tol),
      break;
  end; % checking convergence
end
psi=psi_best;
%figure,surf(psierrM)
if sum(isnan(psi(:))) ~= 0
    psi=nanfree_psi;
end
%Step 5: Getting covariance matrix from spectral factors
for k1 = 1:m
    for k2 = 1:m
        gamtmp(k1,k2,:) = real(ifft(squeeze(psi(k1,k2,:))));
    end
end

%Step 6: Getting noise covariance & transfer function (see Example pp. 424)
A0    = gamtmp(:,:,1);
Z     = A0*A0.'; %Noise covariance matrix not multiplied by sampling frequency

H = zeros(m,m,N+1) + 1i*zeros(m,m,N+1);
for k = 1:N+1
    H(:,:,k) = psi(:,:,k)/A0;       %Transfer function
    S(:,:,k) = psi(:,:,k)*psi(:,:,k)'; %Updated cross-spectral density
end

function gp = PlusOperator(g,nchan,nfreq) %#codegen

% This function is for [ ]+operation:
% to take the positive lags & half of the zero lag and reconstitute
% M. Dhamala, UF, August 2006

% g   = transpose(reshape(g, [nchan^2 2*(nfreq-1)]));
% gam = ifft(g);
%
% % taking only the positive lags and half of the zero lag
% gamp  = gam;
% beta0 = 0.5*gam(1,:);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %aux=reshape(triu(reshape(beta0, [nchan nchan])),[1 nchan^2]);
% %gamp(1,          :) = aux(1,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamp(1,          :) = reshape(triu(reshape(beta0, [nchan nchan])),[1 nchan^2]);
% gamp(nfreq+1:end,:) = 0;
%
% % reconstituting
% gp = fft(gamp);
% gp = reshape(transpose(gp), [nchan nchan 2*(nfreq-1)]);

%------------------------------------------------------
%this is the original code; above is vectorized version
%which is assumed to be faster with many channels present
gam=complex(zeros(size(g)));
for k1 = 1:nchan
 for k2 = 1:nchan
   gam(k1,k2,:) = ifft(squeeze(g(k1,k2,:)));
 end
end

% taking only the positive lags and half of the zero lag
gamp  = gam;
beta0 = 0.5*gam(:,:,1); 
gamp(:,:,1) = triu(beta0);  %this is Stau
gamp(:,:,nfreq+1:end) = 0;

% reconstituting
gp=complex(zeros(size(gamp)));
for k1 = 1:nchan
 for k2 = 1:nchan
   gp(k1,k2,:) = fft(squeeze(gamp(k1,k2,:)));
 end
end

function L=my_chol(M) %#codegen

n = length( M );
L = zeros( n, n );
for i=1:n
    L(i, i) = sqrt( M(i, i) - L(i, :)*L(i, :)' );

    for j=(i + 1):n
        L(j, i) = ( M(j, i) - L(i, :)*L(j, :)' )/L(i, i);
    end
end
L=L.';

%    Partially written by João Rodrigues
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