function DTF=DTF(B,A,N,Fs)

if strcmp(B,'np')
    [K1,K2,N] = size(A);
    DTF=zeros(K1,K1,N);
    h=A;
    for k1=1:K1
        DEN=sum(abs(h(k1,:,:)).^2,2);
        for k2=1:K2
            DTF(k1,k2,:)   = (abs(h(k1,k2,:)).^2)./DEN;
            %DTF(k1,k2,:)   = abs(h(k1,k2,:))./sqrt(DEN);not squarred
        end
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
    
    h=zeros(K1,K1,N);
    DTF=zeros(K1,K1,N);

    for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
            atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end

        
        btmp = zeros(K1);
        for k = 1:q+1,
            btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        h(:,:,n)  = atmp\btmp;
    end

    for k1=1:K1
        DEN=sum(abs(h(k1,:,:)).^2,2);
        for k2=1:K2
            DTF(k1,k2,:)   = (abs(h(k1,k2,:)).^2)./DEN;
            %DTF(k1,k2,:)   = abs(h(k1,k2,:))./sqrt(DEN);not squarred
        end
    end
end

%    Written by Jo?o Rodrigues
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