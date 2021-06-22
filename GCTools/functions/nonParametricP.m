function [H Z S freq]=nonParametricP(X,Fs,alt)

seq_length = size(X,2);
[N M] = size(X);                % Length of signal
NFFT=256;
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
                    [FFT1 freq]=my_fft(hamming(finalPos-k+1)'.*X(i,k:finalPos),NFFT,Fs);
                    [FFT2 freq]=my_fft(hamming(finalPos-k+1)'.*X(j,k:finalPos),NFFT,Fs);
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
                    [FFT1 freq]=my_fft(hamming(finalPos-k+1)'.*X(i,k:finalPos),NFFT,Fs);
                    XSii=XSii+FFT1.*conj(FFT1);
                    it=it+1;
                end
            end
            XSii=XSii/it;
            S(i,i,:)=XSii;
        end
    end
end

if exist('codegen')
    try
        funcName=['wilsonSpectralFactor' N '_mex'];
        if exist(funcName)
            eval(['[H, Z, S, psi] = ' funcName '(S,freq);']);
        else
            eval(['codegen wilsonSpectralFactor -args {S,freq} -o ' funcName]);
            eval(['[H, Z, S, psi] = ' funcName '(S,freq);']);
        end
        
    catch exception
        [H, Z, S, psi] = wilsonSpectralFactor(S,freq);
    end
else
    [H, Z, S, psi] = wilsonSpectralFactor(S,freq);
end



function [Y f]=my_fft(y,NFFT,Fs)
L = length(y);
Y = fft(y,NFFT);
Y=Y(1:NFFT/2+1);
f = Fs/2*linspace(0,1,NFFT/2+1);

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