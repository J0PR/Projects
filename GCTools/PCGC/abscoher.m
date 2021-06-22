function coh=abscoher(X,Fs)

%X npoints*nvars
%falta por frequencas
[spect,freqs]=temp_fft(X, Fs);
for i=1:size(spect,2)
    for j=1:size(spect,2)
        spect_i=spect(:,i);
        spect_j=spect(:,j);
        spect_ij=spect_i.*conj(spect_j);
        coh(i,j,:)=abs(spect_ij)./sqrt(abs(spect_i).*abs(spect_j));
    end
end


function [Spect,f]=temp_fft(y, Fs)
T = 1/Fs;
L = size(y,1);
NFFT = 2^nextpow2(L);
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
Spect=2*(Y(1:NFFT/2+1,:));