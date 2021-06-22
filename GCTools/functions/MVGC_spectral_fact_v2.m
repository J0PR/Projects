function [H Z S freq AF]=MVGC_spectral_fact_v2(X,Fs,freqs)
%TESTING
if Fs>0 %using Welch periodogram... this has to be decided or an argument.
    seq_length = size(X,2);
    [N M] = size(X);                % Length of signal
    NFFT=size(freqs,2);
    %ll=floor(M/NFFT);
    %ll=floor(M/ll);
    %ll=Fs*4;
    if Fs<100
        ll=36;
    else
        ll=Fs*2;
    end
    fftlen1=NFFT/2+1;
    %fftlen2=NFFT+1;
    XSij=zeros(1,fftlen1);
    XSji=zeros(1,fftlen1);
    XSii=zeros(1,fftlen1);
    S=zeros(N,N,fftlen1);
    overlap=0.3;
    for i=1:N
        for j=1:i
            if j~=i
                it=0;
                for k=1:floor(ll*overlap):M
                    finalPos=k+ll-1;
                    if finalPos<=M
                        FFT1=DFT(hamming(finalPos-k+1)'.*X(i,k:finalPos),freqs,Fs);
                        FFT2=DFT(hamming(finalPos-k+1)'.*X(j,k:finalPos),freqs,Fs);
                        XSij=XSij+FFT1.*conj(FFT2);
                        XSji=XSji+FFT2.*conj(FFT1);
                        it=it+1;
                    end
                end
                XSij=XSij/it;
                XSji=XSji/it;
                S(i,j,:)=XSij;
                S(j,i,:)=XSji;
            else
                it=0;
                for k=1:floor(ll*overlap):M
                    finalPos=k+ll-1;
                    if finalPos<=M
                        FFT1=DFT(hamming(finalPos-k+1)'.*X(i,k:finalPos),freqs,Fs);
                        XSii=XSii+FFT1.*conj(FFT1);
                        it=it+1;
                    end
                end
                XSii=XSii/it;
                S(i,i,:)=XSii;
            end
        end
    end
else %alternative use multitaper 
    seq_length = size(X,2);
    %time_halfbandwidth = 2.5;
    %num_seq = 2*(2.5)-1;
    num_seq =5;
    time_halfbandwidth=num_seq/2+1;
    %Obtain DPSSs
    [dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);
    tapers=dps_seq';
    [N M] = size(X);                % Length of signal
    %NFFT = 2^nextpow2(M); % Next power of 2 from length of y
    %if NFFT >256
    NFFT=size(freqs,2);
    %end
    XSij=zeros(1,NFFT/2+1);
    XSji=zeros(1,NFFT/2+1);
    XSii=zeros(1,NFFT/2+1);
    S=zeros(N,N,NFFT/2+1);
    for i=1:N
        for j=1:i
            if j~=i
                for k=1:size(tapers,1)
                    FFT1=DFT(tapers(k,:).*X(i,:),freqs,Fs);
                    FFT2=DFT(tapers(k,:).*X(j,:),freqs,Fs);
                    XSij=XSij+FFT1.*conj(FFT2);
                    XSji=XSji+FFT2.*conj(FFT1);
                end
                XSij=XSij/size(tapers,1);
                XSji=XSji/size(tapers,1);
                S(i,j,:)=XSij;
                S(j,i,:)=XSji;
            else
                for k=1:size(tapers,1)
                    FFT1=DFT(tapers(k,:).*X(i,:),freqs,Fs);
                    XSii=XSii+FFT1.*conj(FFT1);
                end
                XSii=XSii/size(tapers,1);
                S(i,i,:)=XSii;
            end
        end
    end
end

warning off all

[H,Z,AF]=spectral_fact(S,2);
H(:,:,end+1)=zeros(size(H,1),size(H,2));


function [H,SIG,AF]=spectral_fact(S,opt)
if nargin ==1
    opt=2;
end
if opt==1
    [H,SIG,iters] = cpsd_to_var(S);
else
    [G,q] = cpsd_to_autocov(S);
    [AF,SIG] = autocov_to_var(G);
    H = var2trfun(AF,size(AF,3));
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


