function [S,freq]=xmtm(X,NFFT,Fs) 

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
    %end
    XSij=zeros(1,NFFT/2+1);
    XSji=zeros(1,NFFT/2+1);
    XSii=zeros(1,NFFT/2+1);
    S=zeros(N,N,NFFT/2+1);
    for i=1:N
        for j=1:i
            if j~=i
                for k=1:size(tapers,1)
                    [FFT1 freq]=local_fft(tapers(k,:).*X(i,:),NFFT,Fs);
                    [FFT2 freq]=local_fft(tapers(k,:).*X(j,:),NFFT,Fs);
                    XSij=XSij+FFT1.*conj(FFT2);
                    XSji=XSji+FFT2.*conj(FFT1);
                end
                XSij=XSij/size(tapers,1);
                XSji=XSji/size(tapers,1);
                S(i,j,:)=XSij;
                S(j,i,:)=XSji;
            else
                for k=1:size(tapers,1)
                    [FFT1 freq]=local_fft(tapers(k,:).*X(i,:),NFFT,Fs);
                    XSii=XSii+FFT1.*conj(FFT1);
                end
                XSii=XSii/size(tapers,1);
                S(i,i,:)=XSii;
            end
        end
    end
    
function [Y f]=local_fft(y,NFFT,Fs)
L = length(y);
Y = fft(y,NFFT);

Y=Y(1:NFFT/2+1);
f = Fs/2*linspace(0,1,NFFT/2+1);