function [GCout, analysis_freqs] = TV_GrangerMetrics(data,optStruct) %incompleto...

%outputs: GCout, analysis_freqs

AR_mode=optStruct.AR_mode;
extended_AR=optStruct.extended_AR;
measure=optStruct.measure;
spect_norm_opt=optStruct.spect_norm_opt;
NLAGS=optStruct.NLAGS;
freqs=optStruct.freqs;
Fs=optStruct.Fs;
spect_opt=optStruct.spect_opt;

analysis_freqs=[];
GCout=[];

[nvars, npoints]=size(data);
fres=0.1;

if strcmp(freqs,'all') || strcmp(freqs,'mean')
    analysis_freqs = [0:fres:Fs/2-fres];
else
    N_freqs=100;
    analysis_freqs=zeros(size(freqs,1),N_freqs);
    for i=1:size(freqs,1)
        analysis_freqs(i,:)=linspace(freqs(i,1),freqs(i,2),N_freqs);
    end
end

switch AR_mode
    case 'VARMA'
        Sstruct=spectrograms(data,spect_norm_opt,Fs);
        [phi,theta,C]=TV_VARMA(data,NLAGS,NLAGS);%TODO: estimate MA coefficients order outside.
        N=size(phi,3);
        M=zeros(nvars,nvars,size(analysis_freqs,2),N);
        for i=1:N
            covM = C(:,:,i);
            if optStruct.VARMA2MVAR % convert VARMA to MVAR with Matlab's Econometrics Toolbox VGXAR function
                phi_cell=coefMtoCell(phi(:,:,i));
                theta_cell=coefMtoCell(theta(:,:,i));
                varma_model = vgxset('AR', phi_cell,'MA', theta_cell);
                mvar_model=vgxar(varma_model,NLAGS);
                ar_cell=mvar_model.AR;
                if size(ar_cell,1)>size(ar_cell,2)
                    ar_cell=ar_cell';
                end
                AR=cell2mat(ar_cell);
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                end
                ar=[eye(size(AR,1)), -AR];%to be used in the following functions when transforming to freq.verificar se ar ja n tem eye.
                ma=eye(size(AR,1));
            else
                AR=phi(:,:,i);
                MA=theta(:,:,i);
                if extended_AR
                    [AR,covM]=extend_AR(AR,covM);
                    [MA,covM]=extend_AR(MA,covM);
                end
                ar=[eye(size(AR,1)), -AR];
                ma=[eye(size(MA,1)), MA];
            end
            m = computeMetric(data,measure,ma,ar,AR,covM,analysis_freqs,freqs,Fs,spect_norm_opt,Sstruct,i);
            if isempty(m)
                res=[];
                return;
            end
            M(:,:,:,i) = m;
        end
    case 'MVAR'
        Sstruct=spectrograms(data,spect_norm_opt,Fs);
        [phi,C]=TV_MVAR(data,NLAGS);
        N=size(phi,3);
        M=zeros(nvars,nvars,size(analysis_freqs,2),N);
        for i=1:N
            AR=phi(:,:,i);
            covM = C(:,:,i);
            if extended_AR
                [AR,covM]=extend_AR(AR,covM);
            end
            ar=[eye(size(AR,1)), -AR];
            ma=eye(size(AR,1));
            m = computeMetric(data,measure,ma,ar,AR,covM,analysis_freqs,freqs,Fs,spect_norm_opt,Sstruct,i);
            if isempty(m)
                res=[];
                return;
            end
            M(:,:,:,i) = m;
        end
    case {'Wilson-Burg', 'autocov-seq'}
        nAvg=optStruct.spect_avg;
        switch spect_opt
            case {'syncWT','WT'}
                [Cx, W, T, w, analysis_freqs] = syncsq_xWT_lin(data, Fs);
                Sstruct.spectrogram=W;
                Sstruct.freqs=analysis_freqs;
            case 'ST'
                [Cx, S, analysis_freqs] = xST(data, Fs);
                Sstruct.spectrogram=S;
                Sstruct.freqs=analysis_freqs;
            otherwise
                error(['spectrum mode ' spect_opt ' does not exist in time-varying.']);
        end
        N=size(Cx,4);
        for i=1:N
            past=i-nAvg+1;
            if past<1
                past=1;
            end
            cx_i=mean(Cx(:,:,:,past:i),4);
            switch AR_mode
                case 'Wilson-Burg'
                    %% Wilson-Burg
                    [H, Z, ~, ~] = wilsonSpectralFactor(cx_i);
                    ar=H; ma='np'; covM=Z;AR=[];
                    %if strcmp(measure,'newGCFreq') || strcmp(measure,'newGCTime')
                    %UNTESTED
                    [AR,NLAGS] = trfun2var(H);
                    AR=reshape(AR,[nvars nvars*NLAGS]);
                    %end
                case 'autocov-seq'
                    %% Autocovariance sequence
                    [G,~] = cpsd_to_autocov(cx_i);
                    [AR,Z] = autocov_to_var(G);
                    H = var2trfun(AR,size(cx_i,3)-1);
                    NLAGS=size(AR,3);
                    ar=H; ma='np'; covM=Z;
            end
            if extended_AR
                [AR,covM]=extend_AR(AR,covM);
                ar = [eye(size(AR,1)),-AR];
                H=TransfMatrix(eye(size(AR,1)),ar,analysis_freqs,Fs);
                ar=H;
            end
            m = computeMetric(data,measure,ma,ar,AR,covM,analysis_freqs,freqs,Fs,spect_norm_opt,Sstruct,i);
            if isempty(m)
                error(['No metric computed.']);
            end
            M(:,:,:,i) = m;
        end
        if strcmp(spect_opt,'syncWT') && strcmp(freqs,'all')
            %% squeezing
            for i=1:size(M,1)
                for j=1:size(M,2)
                    if i~=j
                        m=flipud(squeeze(M(i,j,:,:)));%flipud to invert y so %flipud to invert y so it has increasng scales (to be used in synsq_causal_squeeze_lin and synsq_cwt_squeeze_lin).
                        m0=m(end,:);% 0 freq
                        m=m(1:end-1,:); %remove 0 freq
                        t=[0:N-1]*(1/Fs);
                        if ~isempty(strfind(measure,'PDC')) || ~isempty(strfind(measure,'DTF'))
                            [ms,~,mavg] = synsq_causal_squeeze_lin(m, squeeze(w(j,:,:)), t, 32, []);
                            ms=ms./mavg;
                        else
                            [ms,~,mavg] = synsq_cwt_squeeze_lin(m, squeeze(w(j,:,:)), t, 32, []);
                        end
                        ms=[m0;ms];
                        M(i,j,:,:)=ms;
                    end
                end
            end
            %%
        end
    otherwise
        error(['AR mode ' AR_mode ' does not exist in time-varying.']);
end
GCout=M;

%% function to compute the causal metric and weight it according to spectral content
function GCout = computeMetric(X,measure,ma,ar,AR,covM,analysis_freqs,freqs,Fs,spect_norm_opt,Sstruct,iter)

for i=1:size(freqs,1)
    switch measure
        case {'DTF', 'PDC', 'DC','oPDC','TransfMatrix'}
            fn=str2func(measure);
            res = fn(ma,ar,analysis_freqs(i,:),Fs);
        case {'dDTF','gDTF', 'gPDC','goPDC'}
            fn=str2func(measure);
            res = fn(ma,ar,covM,analysis_freqs(i,:),Fs);
        case {'newGCFreq'}
            for j=1:size(X,1)
                Sj=abs(stran(X(j,:)));
                temp_freqs = linspace(0,Fs/2,size(Sj,1));
                temp_spect(j,:) = interp1(temp_freqs,squeeze(Sj(:,iter)),analysis_freqs);
            end
            res=newGCFreq(X,AR,covM,analysis_freqs(i,:),Fs,temp_spect);
        case {'pmGGC'}
            res = pmGGC(ma,ar,covM,analysis_freqs(i,:),Fs);
        otherwise
            error(['Function ' measure ' does not exist.']);
    end
    if isempty(Sstruct.spectrogram)
        spectrum=[];
    else
        spectrum=Sstruct.spectrogram(:,:,iter);
    end
    
    res = spectralWeighting(res,ma,ar,covM,spectrum,spect_norm_opt,Sstruct.freqs,analysis_freqs(i,:),Fs);
    
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
function GCout=spectralWeighting(m,ma,ar,covM,spectrum,spect_norm_opt,freqs,analysis_freqs,Fs)

if isempty(spectrum) || isempty(spect_norm_opt) || strcmp(spect_norm_opt,'')
    GCout=m;
    return;
end

switch spect_norm_opt
    case {'TransfMatrix'} %Transfer matrix
        H = TransfMatrix(ma,ar,analysis_freqs,Fs);
        GCout=m.*H;
        return;
    case {'ARSpectrum'} %AR spectral estimate
        spectrum=abs(ARSpect(ma,ar,covM,analysis_freqs,Fs));
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
    case 'none'
        GCout=m;
        return;
    otherwise
        if isempty(spectrum)
            GCout=m;
            return;
        end
end
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



%% pre spectral wighting calculations
function Spectrogram=spectrograms(data,spect_norm_opt,Fs)
if isempty(spect_norm_opt)
    Spectrogram.spectrogram=[];
    Spectrogram.freqs=[];
    return;
end
switch spect_norm_opt
    case {'HHT'} %Hilbert Huang Spectrum
        HTopts.norm_opt = 'NHT';
        HTopts.Fs=Fs;
        HTopts.Fres=0.1;
        
        EMDopts.samplesOut= 0;
        EMDopts.emd_func='NA-EMD';
        EMDopts.N=100;
        EMDopts.Nstd=26;
        EMDopts.stat_val=1;
        EMDopts.validate_set_flag=1;
        EMDopts.pval=0.05;
        [spectrogram,freqs]=data2HHTSpectrogram(data,EMDopts,HTopts);
        
    case {'syncWT'} %Synchrosqueezed Wavelet Transform
        t=[0:size(data,2)-1]*(1/Fs);
        for i=1:size(data,1)
            [spectrogram(i,:,:), freqs, ~, ~, ~]=synsq_cwt_fw_lin(t,data(i,:));
        end
    case {'WT'} %Wavelet Transform
        t=[0:size(data,2)-1]*(1/Fs);
        for i=1:size(data,1)
            [~, freqs, spectrogram(i,:,:), ~, ~]=synsq_cwt_fw_lin(t,data(i,:));
        end
    case {'ST'} % Stockwell Transform
        for i=1:size(data,1)
            spectrogram(i,:,:)=stran(data(i,:));
        end
        freqs = linspace(0,Fs/2,size(Si,1));
    otherwise
        Spectrogram.spectrogram=[];
        Spectrogram.freqs=[];
        return;
end
Spectrogram.spectrogram=abs(spectrogram);
Spectrogram.freqs=freqs;


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