function GGC=NTrialsGewekeGC(X,B,A,C,NLAGS,N,Fs,Nr,Nl)
%LSGewekeGC so sub-model is computed with OLS AR model...
if strcmp(B,'np')
    [K1,K2,N] = size(A);
    h=A;
    if K1==2
        S=NLAGS;
        [P_l P_u]=calculatePMatrix(C);
        corr_cov_l=P_l*C*P_l.'.*eye(size(P_l));
        corr_cov_u=P_u*C*P_u.'.*eye(size(P_u));
        
        corr_h_l=rightDivision(h,P_l);
        corr_h_u=rightDivision(h,P_u);
        for n=1:size(S,3)
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
        %atmp=AL(:,:,n);
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
    GGC=conditionalGC(h,C,X,Fs,NLAGS,N,K2,K1,B,Nr,Nl);
else
    for k1=1:K1;
        for k2=1:K2;
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

function GC=conditionalGC(condH,condC,data,Fs,NLAGS,N,from,to,np_FLAG,Nr,Nl)

dim_xyz = [1, 1, from-2];


GC=NaN*ones(to,from,N);
for i=1:from
    X=data([1:i-1 i+1:end],:);
    if strcmp(np_FLAG,'np')
        [h C S]=calcHandCovNParam(X,Fs,1);
        [K1,K2] = size(h);
    else
        [AR,PE]=armorf(X,Nr,Nl,NLAGS);
        M = size(AR,1);
        A = [eye(M),-AR]; B = eye(M); C  = PE;
        
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
        
        %%%%%%
%         Ashaped=reshape(A(:,K1+1:end),K1,K1,p);
%         AL = calcAfHf(Ashaped, N, Fs, f);
        %h=H;
        %%%%%%
        
        for n=1:N,
            atmp = zeros(K1);
            for k = 1:p+1,
                atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
            end;
%             atmp=AL(:,:,n);
            btmp = zeros(K1);
            for k = 1:q+1,
                btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
            end;
            h(:,:,n)  = atmp\btmp;
        end
    end
    for j=1:to
        if j~=i
            [tempCondH tempCondC]=HInd(condH, condC, i, j);
            [tempH tempC]=subHInd(h, C, i, j);
            %Pcond=calculatePMatrix(tempCondC); ver s é igual a de baixo
            Pcond=calcP_lower(tempCondC,dim_xyz);
            %Psub=calculatePMatrix(tempC); é feito ah frente
            
            tempCondC=Pcond*tempCondC*Pcond.';
            %tempC=Psub*tempC*Psub.';
            
            tempCondH=rightDivision(tempCondH,Pcond);
            %tempH=rightDivision(tempH,Psub); é feito ah frente
            
            S_xx_bar = zeros(dim_xyz(1),dim_xyz(1),N);
            f_ytox_given_z = zeros(1,N);
            G = zeros(dim_xyz(1)+dim_xyz(3),dim_xyz(1)+dim_xyz(3),N);
            Q = zeros(sum(dim_xyz),sum(dim_xyz),N);
            
            Psub = [eye(dim_xyz(1)) zeros(dim_xyz(1),dim_xyz(3));...
            -tempC(1:dim_xyz(1),dim_xyz(1)+1:end)'/(tempC(1:dim_xyz(1),1:dim_xyz(1))) eye(dim_xyz(3))];
        
            for ii = 1:N
                
                G(:,:,ii) = tempH(:,:,ii)/(Psub);
                Q(:,:,ii) = ([G(1:dim_xyz(1),1:dim_xyz(1),ii),zeros(dim_xyz(1),dim_xyz(2)),G(1:dim_xyz(1),dim_xyz(1)+1:end,ii);...
                    zeros(dim_xyz(2),dim_xyz(1)),   eye(dim_xyz(2)),  zeros(dim_xyz(2),dim_xyz(3));...
                    G(dim_xyz(1)+1:end,1:dim_xyz(1),ii), zeros(dim_xyz(3),dim_xyz(2)), G(dim_xyz(1)+1:end,dim_xyz(1)+1:end,ii)])\tempCondH(:,:,ii);
                S_xx_bar(:,:,ii) = Q(1:dim_xyz(1),1:dim_xyz(1),ii)*tempCondC(1:dim_xyz(1),1:dim_xyz(1))*Q(1:dim_xyz(1),1:dim_xyz(1),ii)'+...
                    Q(1:dim_xyz(1),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),ii)*tempCondC(dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2))*Q(1:dim_xyz(1),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),ii)'+...
                    Q(1:dim_xyz(1),dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),ii)*tempCondC(dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz))*Q(1:dim_xyz(1),dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),ii)';
                f_ytox_given_z(ii) = log(det(S_xx_bar(:,:,ii))/det(Q(1:dim_xyz(1),1:dim_xyz(1),ii)*tempCondC(1:dim_xyz(1),1:dim_xyz(1))*Q(1:dim_xyz(1),1:dim_xyz(1),ii)'));
            end
            
            
            %             Q=calcQ(tempH,tempCondH,2);
            %             S=zeros(to,from,N);
            %             for n=1:N
            %                 %S(:,:,n)=Q(:,:,n)*condC*Q(:,:,n)'/Fs;
            %                 S(:,:,n)=Q(:,:,n)*tempCondC*Q(:,:,n)';
            %             end
            %             GC(j,i,:)=log(abs(S(1,1,:))./(tempCondC(1,1)*abs(Q(1,1,:)).^2));
            
            GC(j,i,:)=f_ytox_given_z;
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
cov(:,1:N)=cov(:,arrayAux);
cov(1:N,:)=cov(arrayAux,:);
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