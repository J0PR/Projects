function fGC=newGCFreq(X,a,covM,freqs,Fs,opt_Xf)
% X (nvars*N)
% a (nvars*(nvars*m))
% covM (nvars*nvars)
% freqs (1*nfreqs)
% Fs 1*1

% This method follows the square root of Eq. 30 from [1] or Eq. 4 from [2].
% In Eq. 4 from [2] the ratio is the square root of Eq. 30 in [1]. Ref [2]
% is chosen as it is more recent (2012 vs. 2011, among other reasons..)
% Refs:
% [1] S. Hu et al., "Causality analysis of neural connectivity: critical examination of
% existing methods and advances of new methods.", IEEE transactions on neural networks, 2011
%
% [2] S. Hu, H. Liang, "Causality analysis of neural connectivity: New tool
% and limitations of spectral Granger causality", Neurocomputing, 2012

[nvars N]=size(X);
[nvars nvars_m]=size(a);
m=nvars_m/nvars;
%reshape (nvars*(nvars*P))->((nvars*nvars)*P)
a=reshape(a,[nvars nvars m]);

NFFT = length(freqs);
z = -1i*2*pi/Fs;

fGC=zeros(nvars,nvars,NFFT);
A=zeros(nvars,nvars,NFFT);
spectX=zeros(nvars,NFFT);
if nargin < 6 % Xf is not sent as argument
    Xf=zeros(nvars,NFFT);
else
    Xf=opt_Xf;
end
for n=1:NFFT
        atmp = zeros(nvars);
        for k = 1:m
            atmp = atmp + a(:,:,k)*exp(z*freqs(n)*k);
        end;
        A(:,:,n)  = atmp;
        if nargin < 6 % Xf is not sent as argument
            xtmp=zeros(nvars,1);
            for k = 1:N
                xtmp = xtmp + X(:,k)*exp(z*freqs(n)*k);
            end
            %Xf(:,n) = abs(xtmp/N).^2;
            Xf(:,n) = xtmp/N;
        end
        spectX(:,n)=abs(Xf(:,n)).^2;
        for i=1:nvars
            for k=1:nvars
                if i~=k
                    numerator = abs(A(k,i,n))^2*spectX(i,n);
  
                    den1=0;
                    den2=covM(k,k);
                    %residual array was replaced by covariance matrix.
                    %den2=var(er(k,:));
                    for h=1:nvars
                        if h~=k
                            den1=den1+abs(A(k,h,n))^2*spectX(h,n);
                        end
                    end
                    fGC(k,i,n)=numerator/(den1+den2); %Eq. 30 from [1]. 
                    %fGC(k,i,n)=sqrt(numerator/(den1+den2)); %Do sqrt and get Eq. 4 from [2].
                end
            end
        end
end

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
