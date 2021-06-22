function GGC=pmGGC(B,A,C,N,Fs)

if strcmp(B,'np')
    [K1,K2,N] = size(A);
    h=A;
    if K1==2
        [P_l P_u]=calculatePMatrix(C);
        corr_cov_l=P_l*C*P_l.'.*eye(size(P_l));
        corr_cov_u=P_u*C*P_u.'.*eye(size(P_u));
        
        corr_h_l=rightDivision(h,P_l);
        corr_h_u=rightDivision(h,P_u);
        for n=1:size(h,3)
            S_l(:,:,n)  = corr_h_l(:,:,n)*corr_cov_l*corr_h_l(:,:,n)';
            S_u(:,:,n)  = corr_h_u(:,:,n)*corr_cov_u*corr_h_u(:,:,n)';
        end
    end
else
    [K1,K2] = size(A);
    p = K2/K1-1;
    %a=ones(1,p+1);
    [K1,K2] = size(B);
    q = K2/K1-1;
    %b=ones(1,q+1);
    
    if all(size(N)==1),
        f = (0:N-1)*(Fs/(2*N));
    else
        f = N;
        N = length(N);
    end;
    z = 1i*2*pi/Fs;
    
    h=zeros(K1,K1,N);
    corr_h_l=zeros(K1,K1,N);
    corr_h_u=zeros(K1,K1,N);
    S=zeros(K1,K1,N);
    GGC=zeros(K1,K1,N);
    
    %normalização segundo Geweke
    [P_l P_u]=calculatePMatrix(C);
    corr_cov_l=P_l*C*P_l.'.*eye(size(P_l));
    corr_cov_u=P_u*C*P_u.'.*eye(size(P_u));
    
     %%%%%%
%     Ashaped=reshape(A(:,K1+1:end),K1,K1,p);
%     AL = calcAfHf(Ashaped, N, Fs, f);
    %h=H;
    %%%%%%
    
    for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
            atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        %%%%%%
%         atmp=AL(:,:,n);
        %%%%%%
        corr_atmp_l=P_l*atmp;
        corr_atmp_u=P_u*atmp;
        
        btmp = zeros(K1);
        for k = 1:q+1,
            btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        corr_h_l(:,:,n)  = corr_atmp_l\btmp;
        corr_h_u(:,:,n)  = corr_atmp_u\btmp;
        h(:,:,n)  = atmp\btmp;
        %S(:,:,n)  = h(:,:,n)*C*h(:,:,n)'/Fs;
        S(:,:,n)  = h(:,:,n)*C*h(:,:,n)';
        S_l(:,:,n)  = corr_h_l(:,:,n)*corr_cov_l*corr_h_l(:,:,n)';
        S_u(:,:,n)  = corr_h_u(:,:,n)*corr_cov_u*corr_h_u(:,:,n)';
    end
end

if K1>2
    GGC=conditionalGCPartitions(h,C,N,K1,K1);
else
    for k1=1:K1;
        for k2=1:K1;
            if k2>k1 %lower triang
                GGC(k1,k2,:)=log(abs(S_l(k1,k1,:))./(corr_cov_l(k1,k1)*abs(corr_h_l(k1,k1,:)).^2));
            elseif k2<k1 %upper triang
                GGC(k1,k2,:)=log(abs(S_u(k1,k1,:))./(corr_cov_u(k1,k1)*abs(corr_h_u(k1,k1,:)).^2));
            end
        end
    end
end
end

function Q=calcQ(G,H,index)
res=zeros(size(H));
Q=zeros(size(H));
res(1:index-1,1:index-1,:)=G(1:index-1,1:index-1,:);
res(index+1:end,1:index-1,:)=G(index:end,1:index-1,:);
res(1:index-1,index+1:end,:)=G(1:index-1,index:end,:);
res(index+1:end,index+1:end,:)=G(index:end,index:end,:);
res(index,index,:)=1;
for i=1:size(res,3)
    Q(:,:,i)=res(:,:,i)\H(:,:,i);
end
end

function GC=conditionalGCPartitions(condH,condC,N,from,to)

model_dim = [1, 1, from-2];

GC=NaN*ones(to,from,N);
for i=1:from
    K1=size(condH,1)-1;
    %calculo da matriz de covariancia da particao da matriz H
    %[pC Hpart]=partitionCov(condH,condC,i);
    for j=1:to
        if j~=i
            [tempCondH tempCondC]=HInd(condH, condC, i, j);
            [tempH tempC]=partitionImprov(tempCondH, tempCondC,model_dim);
            %[tempH tempC]=subHInd(Hpart, pC, i, j);
            %Pcond=calculatePMatrix(tempCondC);
            %Psub=calculatePMatrix(tempC); done later
            Pcond=calcP_lower(tempCondC,model_dim);
            tempCondC=Pcond*tempCondC*Pcond.';
            %tempC=Psub*tempC*Psub.'.*eye(size(Psub));
            
            tempCondH=rightDivision(tempCondH,Pcond);
            
            Sxx = zeros(model_dim(1),model_dim(1),N);
            F = zeros(1,N);
            G = zeros(model_dim(1)+model_dim(3),model_dim(1)+model_dim(3),N);
            Q = zeros(sum(model_dim),sum(model_dim),N);
            
            for ii = 1:N
                Psub = [eye(model_dim(1)) zeros(model_dim(1),model_dim(3));...
             -tempC(1:model_dim(1),model_dim(1)+1:end,ii)'/(tempC(1:model_dim(1),1:model_dim(1),ii)) eye(model_dim(3))];
         G(:,:,ii) = tempH(:,:,ii)/(Psub);
         Q(:,:,ii) = ([G(1:model_dim(1),1:model_dim(1),ii),zeros(model_dim(1),model_dim(2)),G(1:model_dim(1),model_dim(1)+1:end,ii);...
             zeros(model_dim(2),model_dim(1)),   eye(model_dim(2)),  zeros(model_dim(2),model_dim(3));...
             G(model_dim(1)+1:end,1:model_dim(1),ii), zeros(model_dim(3),model_dim(2)), G(model_dim(1)+1:end,model_dim(1)+1:end,ii)])\tempCondH(:,:,ii);
         Sxx(:,:,ii) = Q(1:model_dim(1),1:model_dim(1),ii)*tempCondC(1:model_dim(1),1:model_dim(1))*Q(1:model_dim(1),1:model_dim(1),ii)'+...
             Q(1:model_dim(1),model_dim(1)+1:model_dim(1)+model_dim(2),ii)*tempCondC(model_dim(1)+1:model_dim(1)+model_dim(2),model_dim(1)+1:model_dim(1)+model_dim(2))*Q(1:model_dim(1),model_dim(1)+1:model_dim(1)+model_dim(2),ii)'+...
             Q(1:model_dim(1),model_dim(1)+model_dim(2)+1:sum(model_dim),ii)*tempCondC(model_dim(1)+model_dim(2)+1:sum(model_dim),model_dim(1)+model_dim(2)+1:sum(model_dim))*Q(1:model_dim(1),model_dim(1)+model_dim(2)+1:sum(model_dim),ii)';
         F(ii) = log(det(Sxx(:,:,ii))/det(Q(1:model_dim(1),1:model_dim(1),ii)*tempCondC(1:model_dim(1),1:model_dim(1))*Q(1:model_dim(1),1:model_dim(1),ii)'));
            end
            
            %             for jj=1:size(Psub,3)
            %                 tempH(:,:,jj)=tempH(:,:,jj)/Psub(:,:,jj);
            %             end
            
            %             Q=calcQ(tempH,tempCondH,2);
            %             S=zeros(to,from,N);
            %             for n=1:N
            %                 %S(:,:,n)=Q(:,:,n)*condC*Q(:,:,n)'/Fs;
            %                 S(:,:,n)=Q(:,:,n)*tempCondC*Q(:,:,n)';
            % %             end
            %             GC(j,i,:)=log(abs(S(1,1,:))./(tempCondC(1,1)*abs(Q(1,1,:)).^2));
            GC(j,i,:)=abs(F);
        end
    end
    
end
end

function [H cov]=HInd(H, cov, from, to)
N=size(H,1);
remaining=find(~ismember(1:N, [to from]));
arrayAux=[to from remaining];

cov(:,1:N)=cov(:,arrayAux);
cov(1:N,:)=cov(arrayAux,:);
for i=1:size(H,3)
    H(:,1:N,i)=H(:,arrayAux,i);
    H(1:N,:,i)=H(arrayAux,:,i);
end
end

function [H cov]=subHInd(H, cov, from, to)
N=size(H,1);
if to>from
    to=to-1;
end
remaining=find(~ismember(1:N, [to]));
arrayAux=[to remaining];
%arrayAux=1:N;
%arrayAux(1)=to;
%arrayAux(to)=1;
for i=1:size(cov,3)
    cov(:,1:N,i)=cov(:,arrayAux,i);
    cov(1:N,:,i)=cov(arrayAux,:,i);
end
for i=1:size(H,3)
    H(:,1:N,i)=H(:,arrayAux,i);
    H(1:N,:,i)=H(arrayAux,:,i);
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