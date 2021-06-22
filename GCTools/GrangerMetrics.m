function [GCout, analysis_freqs] = GrangerMetrics(X,optStruct)

%outputs: GCout, analysis_freqs

AR_mode=optStruct.AR_mode;
if ~isfield(optStruct,'extended_AR')
    optStruct.extended_AR=0;
end
if ~isfield(optStruct,'spect_opt')
    optStruct.spect_opt='';
end
extended_AR=optStruct.extended_AR;
measure=optStruct.measure;
NLAGS=optStruct.NLAGS;
freqs=optStruct.freqs;
Fs=optStruct.Fs;
PVAL=optStruct.pval;
spect_opt=optStruct.spect_opt;

analysis_freqs=[];
GCout=[];

%[nvars, npoints]=size(X);
%fres=0.1;

%% GCI
if strcmp(measure, 'GCI')
    switch AR_mode
        case {'MVAR_Burg'}
            ret = tGC_mvar(X,NLAGS,extended_AR);
            [PR,q] = bonferroni_significance(ret,PVAL);
            if PVAL ~=0 && PVAL ~=1
                GCout = ret.gc.*PR;
            else
                GCout = ret.gc;
            end
        case {'MVAR_OLS'}
            ret = cca_granger_regress(X,NLAGS,1);
            [PR,q] = bonferroni_significance(ret,PVAL);
            if PVAL ~=0 && PVAL ~=1
                GCout = ret.gc.*PR;
            else
                GCout = ret.gc;
            end
        case {'eMVAR_OLS'}
            ret = tGC_0lag_out(X,NLAGS);
            [PR,q] = bonferroni_significance(ret,PVAL);
            if PVAL ~=0 && PVAL ~=1
                GCout = ret.gc.*PR;
            else
                GCout = ret.gc;
            end
        otherwise
            error(['AR mode ' AR_mode ' does not exist for GCI.']);
    end
    return;
end

if strcmp(freqs,'all') || strcmp(freqs,'mean')
    %analysis_freqs = [0:fres:Fs/2-fres];
    N_freqs=100;
    analysis_freqs = linspace(0,Fs/2,N_freqs);
else
    N_freqs=100;
    analysis_freqs=zeros(size(freqs,1),N_freqs);
    for i=1:size(freqs,1)
        analysis_freqs(i,:)=linspace(freqs(i,1),freqs(i,2),N_freqs);
    end
end

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
        AR=reshape(phi,[size(theta,1) size(theta,2)*size(theta,3)]);
        MA=reshape(theta,[size(theta,1) size(theta,2)*size(theta,3)]);
        
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
            %ma = eye(size(AR,1));
        end
        %er=strv.lastResid;
    case {'Wilson-Burg'} %REMOVE WHEN V2 READY
        %% Wilson-Burg
        [Cx,analysis_freqs,w]=xSpect(X,Fs,optStruct);
        [H, Z, ~, ~] = wilsonSpectralFactor(Cx);
        ar=H; ma='np'; covM=Z;
        %if strcmp(measure,'newGCFreq') || strcmp(measure,'newGCTime')
        %UNTESTED
        [AR,NLAGS] = trfun2var(H);
        M = size(AR,1);
        AR=reshape(AR,[M M*NLAGS]);
        if extended_AR
            [AR,covM]=extend_AR(AR,covM);
            ar = [eye(M),-AR];
            H=TransfMatrix(eye(M),ar,analysis_freqs,Fs);
            ar=H;
        end
        %[~,E] = var_to_tsdata(AR,Z,npoints);
        %er=E';
        %end
    case {'autocov-seq'} %REMOVE WHEN V2 READY
        %% Autocovariance sequence
        [Cx,analysis_freqs,w]=xSpect(X,Fs,optStruct);
        %analysis_freqs(1)=[];
        [G,~] = cpsd_to_autocov(Cx);
        [AR,Z] = autocov_to_var(G);
        H = var2trfun(AR,size(Cx,3)-1);
        NLAGS=size(AR,3);
        M = size(AR,1);
        ar=H; ma='np'; covM=Z;
        if extended_AR
            [AR,covM]=extend_AR(AR,covM);
            ar = [eye(M),-AR];
            H=TransfMatrix(eye(M),ar,analysis_freqs,Fs);
            ar=H;
        end
    case {'NTrials_MVAR_Burg'}
        %% ARMORF
        [AR,PE]=armorf(X,optStruct.Nr,optStruct.Nl,NLAGS);
        AR=-AR;
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
    
%% newGCTime
if strcmp(measure, 'newGCTime')
    GCout=newGCTime(X,AR,covM);
    return;
end
%%

for i=1:size(freqs,1)
%%% ACTIVATE WHEN V2 IS READY.
%   NP analysis computed for each freq band in analysis_freqs
%     S=[];
%     switch AR_mode
%         case {'Wilson-Burg'}
%             %% Wilson-Burg
%             [H, Z, S]=calcHandCovNParam_v2(X,Fs,analysis_freqs(i,:));
%             ar=H; ma='np'; covM=Z;
%             if strcmp(measure,'newGCFreq')
%                 %UNTESTED
%                 [AR,NLAGS] = trfun2var(H);
%                 %[~,E] = var_to_tsdata(AR,Z,npoints);
%                 %er=E'; activate if needed
%             end
%         case {'MVGC'}
%             %% Autocovariance sequence
%             [H, Z, S, AR]=MVGC_spectral_fact_v2(X,Fs,analysis_freqs(i,:));
%             NLAGS=size(AR,3);
%             ar=H; ma='np'; covM=Z;
%             if strcmp(measure,'newGCFreq')
%                 %UNTESTED
%                 %er=predError(X,AR,NLAGS)'; activate if needed
%             end
%         otherwise
%             error(['Cannot fit data with ' AR_mode '. Method does not exist.']);
%             return;
%     end
    %% Metrics
    switch measure
        case {'DTF', 'ffDTF', 'PDC', 'DC','oPDC','TransfMatrix','rPDC','roPDC'}
            fn=str2func(measure);
            res = fn(ma,ar,analysis_freqs(i,:),Fs);
        case {'dDTF','gDTF', 'gPDC','goPDC','rgPDC','rgoPDC'}
            fn=str2func(measure);
            res = fn(ma,ar,covM,analysis_freqs(i,:),Fs);
        case {'newGCFreq'}
            res=newGCFreq(X,AR,covM,analysis_freqs(i,:),Fs);
        case {'pmGGC','GGC'}
            if strcmp(measure,'GGC')
                res = GGC(X,ma,ar,covM,NLAGS,analysis_freqs(i,:),Fs,optStruct);%need to send optStruct due to the estimation of reduced AR models inside
            else
                res = pmGGC(ma,ar,covM,analysis_freqs(i,:),Fs);
            end
        otherwise
            error(['Function ' measure ' does not exist.']);
    end
    res = spectralWeighting(res,ma,ar,covM,X,optStruct,analysis_freqs(i,:),Fs);
    if strcmp(freqs,'all')
        GCout=res;
        return;
    else
        if strcmp(AR_mode,'Wilson-Burg') || strcmp(AR_mode,'autocov-seq')
            GCout(:,:,i)=mean(res(:,:,analysis_freqs>=freqs(i,1) & analysis_freqs<freqs(i,2)),3);
        else
            GCout(:,:,i)=mean(res,3);
        end
    end
end

%% Spectral Weighting
function GCout=spectralWeighting(m,ma,ar,covM,data,optStruct,analysis_freqs,Fs)
if ~isfield(optStruct,'spect_norm_opt')
    GCout=m;
    return;
elseif isempty(optStruct.spect_norm_opt)
    GCout=m;
    return;
else
    spect_norm_opt=optStruct.spect_norm_opt;
end

switch spect_norm_opt
    case {'FFT'} %FFT spectrum
        if isfield(optStruct,'Nr')
            pos=1;
            spectrum=[];
            max_Nl=max(optStruct.Nl);
            for it=1:optStruct.Nr
                if length(optStruct.Nl)>1
                    temp_len=optStruct.Nl(it);
                else
                    temp_len=optStruct.Nl;
                end
                if pos+max_Nl-1<=size(data,2)
                    tspan=pos:pos+max_Nl-1;
                    temp_data=data(:,tspan);
                else
                    tspan=pos:size(data,2);
                    temp_data=[data(:,tspan) zeros(size(data,1),max_Nl-length(tspan))];
                end
                
                pos=pos+temp_len;
                temp_spectrum=abs(DFT(temp_data,analysis_freqs,Fs));
                if isempty(spectrum)
                    spectrum=temp_spectrum;
                else
                    spectrum=spectrum+temp_spectrum;
                end
            end
            spectrum=spectrum/optStruct.Nr;
        else
            spectrum=abs(DFT(data,analysis_freqs,Fs));
        end
        spectrum=spectrum./repmat(max(spectrum,[],2),[1 size(spectrum,2)]);
        for j=1:size(m,1)
            spect_j = spectrum(j,:);
            for i=1:size(m,2)
                if i~=j
                    m(i,j,:)=squeeze(m(i,j,:)).*spect_j';
                else
                    m(i,j,:)=spect_j;
                end
            end
        end
        GCout=m;
        return;
    case {'TransfMatrix'} %Transfer matrix
        H = TransfMatrix(ma,ar,analysis_freqs,Fs);
        GCout=m.*H;
        return;
    case {'ARSpectrum'} %AR spectral estimate
        spectrum=ARSpect(ma,ar,covM,analysis_freqs,Fs);
        spectrum=abs(spectrum);
        %spectrum=spectrum/max(max(spectrum));
        spectrum=spectrum./repmat(max(spectrum,[],2),[1 size(spectrum,2)]);
        for j=1:size(m,1)
            spect_j = spectrum(j,:);
            for i=1:size(m,2)
                if i~=j
                    m(i,j,:)=squeeze(m(i,j,:)).*spect_j';
                else
                    m(i,j,:)=spect_j;
                end
            end
        end
        GCout=m;
        return;
    case {'HHT'} %Hilbert Huang Spectrum
        HTopts.norm_opt = 'NHT';
        HTopts.Fs=Fs;
        HTopts.Fres=0.1;
        
        EMDopts.samplesOut= 10;
        EMDopts.emd_func='NA-EMD';
        EMDopts.N=100;
        EMDopts.Nstd=26;
        EMDopts.stat_val=1;
        EMDopts.validate_set_flag=1;
        EMDopts.pval=0.05;
        [spectrum,freqs]=data2HHTSpectrum(data,EMDopts,HTopts);

    case {'syncWT'} %Synchrosqueezed Wavelet Transform
        t=[0:size(data,2)-1]*(1/Fs);
        for i=1:size(data,1)
            [Ti, freqs, ~, ~, ~]=synsq_cwt_fw_lin(t,data(i,:));
            spectrum(i,:)=mean(abs(Ti),2);
        end
    case {'WT'} %Wavelet Transform
        if isfield(optStruct,'Nr')
            pos=1;
            spectrum=[];
            max_Nl=max(optStruct.Nl);
            for it=1:optStruct.Nr
                if length(optStruct.Nl)>1
                    temp_len=optStruct.Nl(it);
                else
                    temp_len=optStruct.Nl;
                end
                if pos+max_Nl-1<=size(data,2)
                    tspan=pos:pos+max_Nl-1;
                    temp_data=data(:,tspan);
                else
                    tspan=pos:size(data,2);
                    temp_data=[data(:,tspan) zeros(size(data,1),max_Nl-length(tspan))];
                end
                
                pos=pos+temp_len;
                %%
                t=[0:size(temp_data,2)-1]*(1/Fs);
                for i=1:size(data,1)
                    [~, freqs, Wi, ~, ~]=synsq_cwt_fw_lin(t,temp_data(i,:));
                    temp_spectrum(i,:)=mean(abs(Wi(:,1:temp_len)),2);
                end
                %%
                if isempty(spectrum)
                    spectrum=temp_spectrum;
                else
                    spectrum=spectrum+temp_spectrum;
                end
                clear temp_spectrum;
            end
            spectrum=spectrum/optStruct.Nr;
        else
            t=[0:size(data,2)-1]*(1/Fs);
            for i=1:size(data,1)
                [~, freqs, Wi, ~, ~]=synsq_cwt_fw_lin(t,data(i,:));
                spectrum(i,:)=mean(abs(Wi),2);
            end
        end
    case {'ST'} % Stockwell Transform
        if isfield(optStruct,'Nr')
            pos=1;
            spectrum=[];
            max_Nl=max(optStruct.Nl);
            for it=1:optStruct.Nr
                if length(optStruct.Nl)>1
                    temp_len=optStruct.Nl(it);
                else
                    temp_len=optStruct.Nl;
                end
                if pos+max_Nl-1<=size(data,2)
                    tspan=pos:pos+max_Nl-1;
                    temp_data=data(:,tspan);
                else
                    tspan=pos:size(data,2);
                    temp_data=[data(:,tspan) zeros(size(data,1),max_Nl-length(tspan))];
                end
                
                pos=pos+temp_len;
                %%
                for i=1:size(temp_data,1)
                    [Si]=stran(temp_data(i,:));
                    temp_spectrum(i,:)=mean(abs(Si(:,1:temp_len)),2);
                end
                %%
                if isempty(spectrum)
                    spectrum=temp_spectrum;
                else
                    spectrum=spectrum+temp_spectrum;
                end
                clear temp_spectrum;
            end
            spectrum=spectrum/optStruct.Nr;
        else
            for i=1:size(data,1)
                [Si]=stran(data(i,:));
                spectrum(i,:)=mean(abs(Si),2);
            end
        end
        
        
        freqs = linspace(0,Fs/2,size(Si,1));
    otherwise
        GCout=m;
        return;
end
spectrum=abs(spectrum);
%spectrum=spectrum/max(max(spectrum));
spectrum=spectrum./repmat(max(spectrum,[],2),[1 size(spectrum,2)]);
for j=1:size(m,1)
    spect_j = interp1(freqs,spectrum(j,:),analysis_freqs);
    for i=1:size(m,2)
        if i~=j
            m(i,j,:)=squeeze(m(i,j,:)).*spect_j';
        else
            m(i,j,:)=spect_j;
        end
    end
end

GCout=m;

%% aux functions
function res=coefMtoCell(coefM)
%coefM: nVars*(nVars*p)
%res: Cell length == p, with nVar*nVar matrices inside
nVars=size(coefM,1);
p=size(coefM,2)/nVars;
if p==1
    res{1}=coefM;
else
    res = mat2cell(coefM,[nVars],repmat(nVars,[1 p]));
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