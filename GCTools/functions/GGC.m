function GGC=GGC(X,B,A,C,NLAGS,N,Fs,optStruct)

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
    corr_h_l=zeros(K1,K1,N);
    corr_h_u=zeros(K1,K1,N);
    S=zeros(K1,K1,N);
    GGC=zeros(K1,K1,N);
    
    [P_l P_u]=calculatePMatrix(C);
    corr_cov_l=P_l*C*P_l.'.*eye(size(P_l));
    corr_cov_u=P_u*C*P_u.'.*eye(size(P_u));
    
    
    for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
            atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;

        corr_atmp_l=P_l*atmp;
        corr_atmp_u=P_u*atmp;
        
        btmp = zeros(K1);
        for k = 1:q+1,
            btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        corr_h_l(:,:,n)  = corr_atmp_l\btmp;
        corr_h_u(:,:,n)  = corr_atmp_u\btmp;
        h(:,:,n)  = atmp\btmp;
        S(:,:,n)  = h(:,:,n)*C*h(:,:,n)';
        S_l(:,:,n)  = corr_h_l(:,:,n)*corr_cov_l*corr_h_l(:,:,n)';
        S_u(:,:,n)  = corr_h_u(:,:,n)*corr_cov_u*corr_h_u(:,:,n)';
    end
end

if K1>2
    GGC=conditionalGC(h,C,X,Fs,NLAGS,N,K1,K1,B,optStruct);
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

function GC=conditionalGC(condH,condC,data,Fs,NLAGS,N,from,to,np_FLAG,optStruct)

AR_mode=optStruct.AR_mode;
extended_AR=optStruct.extended_AR;

model_dim = [1, 1, from-2];
n=from;
GC=NaN*ones(to,from,N);
for i=1:from
    X=data([1:i-1 i+1:end],:);
    if strcmp(np_FLAG,'np')
        switch AR_mode
            case {'Wilson-Burg'} %REMOVE WHEN V2 READY
                %% Wilson-Burg
                [Cx,analysis_freqs,w]=xSpect(X,Fs,optStruct);
                [H, Z, ~, ~] = wilsonSpectralFactor(Cx);
                covM=Z;
                %if strcmp(measure,'newGCFreq') || strcmp(measure,'newGCTime')
                %UNTESTED
                [AR,NLAGS] = trfun2var(H);
                M = size(AR,1);
                AR=reshape(AR,[n n*NLAGS]);
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                    ar = [eye(M),-AR];
                    H=TransfMatrix(eye(M),ar,analysis_freqs,Fs);
                end
                %[~,E] = var_to_tsdata(AR,Z,npoints);
                %er=E';
                %end
            case {'autocov-seq'} %REMOVE WHEN V2 READY
                %% Autocovariance sequence
                [Cx,analysis_freqs,w]=xSpect(X,Fs,optStruct);
                [G,~] = cpsd_to_autocov(Cx);
                [AR,Z] = autocov_to_var(G);
                H = var2trfun(AR,size(AR,3));
                NLAGS=size(AR,3);
                M = size(AR,1);
                covM=Z;
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                    ar = [eye(M),-AR];
                    H=TransfMatrix(eye(M),ar,analysis_freqs,Fs);
                end
                %if strcmp(measure,'newGCFreq')
                %UNTESTED
                %er=predError(X,AR,NLAGS)';
                %end
        end
    else
        
        %% AR MODE selection
        S=[];  %REMOVE WHEN V2 READY or leave it in V3
        %V2: compute spectral fact only for frequencies of interest
        %V3: compute for all frequencies but adjust resolution to have the desired
        %points in the bands(ex: if band with 10Hz width needs 100 points and Nyq freq
        %is 200Hz then the spectrum needs 2000 points.. npoints =
        %(Bpoints/Bwidth)*NyqFreq
        switch AR_mode
            case {'MVAR_Burg'}
                %% MVAR
                [AR,RC,PE,DC,u] = mvar(X',NLAGS);
                M = size(AR,1);
                ar = [eye(M),-AR]; ma = eye(M); covM  = PE(:,(1-M:0)+end);
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                    ar = [eye(M),-AR];
                end
                %er=u';
            case {'MVAR_OLS','eMVAR_OLS'}
                %% LS and CP
                cp_mvar = mvar_0lag_out(X,NLAGS,AR_mode);
                AR=cp_mvar.alpha;
                PE=cp_mvar.Z;
                M = size(AR,1);
                ar = [eye(M),-AR]; ma = eye(M); covM  = PE;
                if extended_AR && ~strcmp(AR_mode,'eMVAR_OLS') %if model already has the AR0, dont do this..
                    [AR,covM]=extend_AR(AR,covM);
                    ar = [eye(M),-AR];
                end
                %er=cp_mvar.u';
            case {'VARMA'}
                %% VARMA
                %MAorder=1;
                strv=VARMA_est(X,NLAGS,NLAGS);
                phi=strv.lastPhi(:,:,2:end); theta=strv.lastTheta(:,:,2:end); covM=strv.lastSigmar;
                AR=reshape(phi,[n n*NLAGS]);
                MA=reshape(theta,[nvars nvars*MAorder]);
                
                if optStruct.VARMA2MVAR % convert VARMA to MVAR with Matlab's Econometrics Toolbox VGXAR function
                    AR_cell=coefMtoCell(AR);
                    MA_cell=coefMtoCell(MA);
                    varma_model = vgxset('AR', AR_cell,'MA', MA_cell);
                    mvar_model=vgxar(varma_model,NLAGS);
                    AR_cell=mvar_model.AR;
                    if size(AR_cell,1)>size(AR_cell,2)
                        AR_cell=AR_cell';
                    end
                    AR=cell2mat(AR_cell);
                    if extended_AR
                        [AR,covM]=extend_AR(AR,covM);
                    end
                    ar=[eye(size(AR,1)) -AR];
                    ma=eye(size(AR,1));
                else
                    if extended_AR
                        [AR,covM]=extend_AR(AR,covM);
                        [MA,covM]=extend_AR(MA,covM);
                    end
                    ar=[eye(size(AR,1)) -AR];
                    ma=[eye(size(AR,1)) MA];
                end
                %er=strv.lastResid;
            case {'NTrials_MVAR_Burg'}
                %% ARMORF
                [AR,PE]=armorf(X,optStruct.Nr,optStruct.Nl,NLAGS);
                M = size(AR,1);
                ar = [eye(M),-AR]; ma = eye(M); covM  = PE;
                %er=predError(X,AR,NLAGS)';
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                    ar = [eye(M),-AR];
                end
            otherwise
                error(['Cannot fit data with ' AR_mode '. Method does not exist.']);
        end
        %%
        
        [AR,RC,PE] = mvar(X',NLAGS);
        M = size(AR,1);
        ar = [eye(M),-AR]; ma = eye(M); covM  = PE(:,(1-M:0)+end);
        
        [K1,K2] = size(ar);
        p = K2/K1-1;
        [K1,K2] = size(ma);
        q = K2/K1-1;
        
        if all(size(N)==1),
            f = (0:N-1)*(Fs/(2*N));
        else
            f = N;
            N = length(N);
        end;
        z = 1i*2*pi/Fs;
        
        H=zeros(K1,K1,N);

        
        for n=1:N,
            atmp = zeros(K1);
            for k = 1:p+1,
                atmp = atmp + ar(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
            end;
            btmp = zeros(K1);
            for k = 1:q+1,
                btmp = btmp + ma(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
            end;
            H(:,:,n)  = atmp\btmp;
        end
    end
    for j=1:to
        if j~=i
            [tempCondH tempCondC]=HInd(condH, condC, i, j);
            [tempH tempC]=subHInd(H, covM, i, j);
            %Pcond=calculatePMatrix(tempCondC); ver s é igual a de baixo
            Pcond=calcP_lower(tempCondC,model_dim);
            %Psub=calculatePMatrix(tempC); é feito ah frente
            
            tempCondC=Pcond*tempCondC*Pcond.';
            %tempC=Psub*tempC*Psub.';
            
            tempCondH=rightDivision(tempCondH,Pcond);
            %tempH=rightDivision(tempH,Psub); é feito ah frente
            
            Sxx = zeros(model_dim(1),model_dim(1),N);
            F = zeros(1,N);
            G = zeros(model_dim(1)+model_dim(3),model_dim(1)+model_dim(3),N);
            Q = zeros(sum(model_dim),sum(model_dim),N);
            
            Psub = [eye(model_dim(1)) zeros(model_dim(1),model_dim(3));...
            -tempC(1:model_dim(1),model_dim(1)+1:end)'/(tempC(1:model_dim(1),1:model_dim(1))) eye(model_dim(3))];
        
            for ii = 1:N
                
                G(:,:,ii) = tempH(:,:,ii)/(Psub);
                Q(:,:,ii) = ([G(1:model_dim(1),1:model_dim(1),ii),zeros(model_dim(1),model_dim(2)),G(1:model_dim(1),model_dim(1)+1:end,ii);...
                    zeros(model_dim(2),model_dim(1)),   eye(model_dim(2)),  zeros(model_dim(2),model_dim(3));...
                    G(model_dim(1)+1:end,1:model_dim(1),ii), zeros(model_dim(3),model_dim(2)), G(model_dim(1)+1:end,model_dim(1)+1:end,ii)])\tempCondH(:,:,ii);
                Sxx(:,:,ii) = Q(1:model_dim(1),1:model_dim(1),ii)*tempCondC(1:model_dim(1),1:model_dim(1))*Q(1:model_dim(1),1:model_dim(1),ii)'+...
                    Q(1:model_dim(1),model_dim(1)+1:model_dim(1)+model_dim(2),ii)*tempCondC(model_dim(1)+1:model_dim(1)+model_dim(2),model_dim(1)+1:model_dim(1)+model_dim(2))*Q(1:model_dim(1),model_dim(1)+1:model_dim(1)+model_dim(2),ii)'+...
                    Q(1:model_dim(1),model_dim(1)+model_dim(2)+1:sum(model_dim),ii)*tempCondC(model_dim(1)+model_dim(2)+1:sum(model_dim),model_dim(1)+model_dim(2)+1:sum(model_dim))*Q(1:model_dim(1),model_dim(1)+model_dim(2)+1:sum(model_dim),ii)';
                F(ii) = log(det(Sxx(:,:,ii))/det(Q(1:model_dim(1),1:model_dim(1),ii)*tempCondC(1:model_dim(1),1:model_dim(1))*Q(1:model_dim(1),1:model_dim(1),ii)'));
            end
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
cov(:,1:N)=cov(:,arrayAux);
cov(1:N,:)=cov(arrayAux,:);
for i=1:size(H,3)
    H(:,1:N,i)=H(:,arrayAux,i);
    H(1:N,:,i)=H(arrayAux,:,i);
end
end

function [Sx,analysis_freqs,w]=xSpect(X,Fs,optStruct)
if isfield(optStruct,'Nr')
    pos=1;
    Sx=[];
    max_Nl=max(optStruct.Nl);
    for it=1:optStruct.Nr
        if length(optStruct.Nl)>1
            temp_len=optStruct.Nl(it);
        else
            temp_len=optStruct.Nl;
        end
        if pos+max_Nl-1<=size(X,2)
            tspan=pos:pos+max_Nl-1;
            temp_X=X(:,tspan);
        else
            tspan=pos:size(X,2);
            temp_X=[X(:,tspan) zeros(size(X,1),max_Nl-length(tspan))];
        end
        pos=pos+temp_len;
        %%
        [temp_Sx,analysis_freqs,w]=xSpect_sub(temp_X,Fs,optStruct.spect_opt,temp_len);
        %%
        if isempty(Sx)
            Sx=temp_Sx;
        else
            Sx=Sx+temp_Sx;
        end
        clear temp_Sx;
    end
    Sx=Sx/optStruct.Nr;
else
    [Sx,analysis_freqs,w]=xSpect_sub(X,Fs,optStruct.spect_opt,size(X,2));
end
end


function [Sx,analysis_freqs,w]=xSpect_sub(X,Fs,spect_opt,temp_len)

Sx=[];
analysis_freqs=[];
w=[];
switch spect_opt
    case {'syncsqWT','WT'}
        [Cx, W, T, w, analysis_freqs] = syncsq_xWT_lin(X, Fs);
        Sx=mean(Cx(:,:,:,1:temp_len),4);
    case 'ST'
        [Cx, St, analysis_freqs] = xST(X, Fs);
        Sx=mean(Cx(:,:,:,1:temp_len),4);
    case 'Welch'
        NFFT=512;
        fftlen=NFFT/2+1;
        Sx=zeros(size(X,1),size(X,1),fftlen);
        for i=1:size(X,1)
            for j=1:i
                [Sx(i,j,:) analysis_freqs]=cpsd(X(i,:),X(j,:),[],[],NFFT,Fs);
                [Sx(j,i,:) analysis_freqs]=cpsd(X(j,:),X(i,:),[],[],NFFT,Fs);
            end
        end
    case 'mtm'
        NFFT=512;
        fftlen=NFFT/2+1;
        [Sx,analysis_freqs]=xmtm(X,NFFT,Fs);
    otherwise
        error(['spectrum mode ' spect_opt ' does not exist in time-varying.']);
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