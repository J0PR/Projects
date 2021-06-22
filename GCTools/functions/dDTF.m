function dDTF=dDTF(B,A,C,N,Fs)
%directed Direct Transfer Function
% Based in Eq. 7 from: 
% A. Korzeniewska,"Determination of information flow direction among brain structures by a
% modified directed transfer function (dDTF) method" 2003

if strcmp(B,'np')
    [K1,K2,N] = size(A);
    dDTF=zeros(K1,K1,N);
    ffDTF=zeros(K1,K1,N);
    pCOH=zeros(K1,K1,N);
    G=zeros(K1,K1,N);
    invC=inv(C);
    h=A;
    for n=1:N
        g=inv(h(:,:,n));
        G(:,:,n) = g'*invC*g;
    end
    for k1=1:K1
        DEN=sum(abs(h(k1,:,:)).^2,2);
        for k2=1:K2
            ffDTF(k1,k2,:) = (abs(h(k1,k2,:)).^2)./sum(DEN,3);
            %ffDTF(k1,k2,:) = abs(h(k1,k2,:))./sqrt(sum(DEN,3));not squarred
            pCOH(k1,k2,:) = abs(G(k1,k2,:).^2)./(G(k1,k1,:).*G(k2,k2,:));
        end
    end
    dDTF = abs(pCOH.*ffDTF); 
else
    [K1,K2] = size(A);
    p = K2/K1-1;
    %a=ones(1,p+1);
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
    dDTF=zeros(K1,K1,N);
    ffDTF=zeros(K1,K1,N);
    pCOH=zeros(K1,K1,N);
    G=zeros(K1,K1,N);
    invC=inv(C);
    
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
        g = atmp/btmp;
        G(:,:,n) = g'*invC*g;
    end
    for k1=1:K1
        DEN=sum(abs(h(k1,:,:)).^2,2);
        for k2=1:K2
            ffDTF(k1,k2,:) = abs((abs(h(k1,k2,:)).^2)./sum(DEN,3));
            %ffDTF(k1,k2,:) = abs(h(k1,k2,:))./sqrt(sum(DEN,3));not squarred
            pCOH(k1,k2,:) = abs(G(k1,k2,:).^2)./(G(k1,k1,:).*G(k2,k2,:));
            
        end
    end
    dDTF = abs(pCOH.*ffDTF);
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