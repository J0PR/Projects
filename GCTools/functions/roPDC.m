function roPDC=roPDC(B,A,N,Fs)
%row-wise orthogonalized PDC
if strcmp(B,'np')
    [K1,K2,N] = size(A);
    roPDC=zeros(K1,K1,N);
    H=A;
    tmp1=zeros(1,K1);
    for n=1:N
        atmp=inv(H(:,:,n));
        for k1 = 1:K1,
            tmp = squeeze(atmp(k1,:))';%row-wise
            tmp1(k1) = sqrt(tmp'*tmp);
        end;
        roPDC(:,:,n)  = (abs(real(atmp)).*abs(imag(atmp)))./(tmp1(ones(1,K1),:)').^2;
    end
else
    [K1,K2] = size(A);
    p = K2/K1-1;
    [K1,K2] = size(B);
    q = K2/K1-1;
    
    if all(size(N)==1),
        f = (0:N-1)*(Fs/(2*N));
    else
        f = N;
        N = length(N);
    end;
    z = 1i*2*pi/Fs;
    
    
    roPDC=zeros(K1,K1,N);
    
    tmp1=zeros(1,K1);
    for n=1:N
        atmp = zeros(K1);
        for k = 1:p+1,
            atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end
        btmp = zeros(K1);
        for k = 1:q+1,
            btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        atmp  = atmp/btmp;
        for k1 = 1:K1,
            tmp = squeeze(atmp(k1,:))';%row-wise
            tmp1(k1) = sqrt(tmp'*tmp);
        end;
        roPDC(:,:,n)  = (abs(real(atmp)).*abs(imag(atmp)))./(tmp1(ones(1,K1),:)').^2;
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