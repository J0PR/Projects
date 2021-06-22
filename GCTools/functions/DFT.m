function Xf=DFT(X,freqs,Fs)
%Discrete Fourier Transform for the desired frequencies
% X: time-series: nvars*npoints
% freqs: array of desired frequencies (from 0 to Nyq freq): 1*numfreqs
% Fs: sampling frequency: 1*1

[nvars, N]=size(X);
z = -1i*2*pi/Fs;
NFFT = length(freqs);
Xf=zeros(nvars,NFFT);
for n=1:NFFT
    xtmp=zeros(nvars,1);
    for k = 1:N
        xtmp = xtmp + X(:,k)*exp(z*freqs(n)*k);
    end
    Xf(:,n) = xtmp/N;
end