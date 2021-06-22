function [AL H] = calcAfHf(A, nf, Fs, f)
%    '''Calculates A(f), in the frequency domain.
%
%     Input:
%         A(nChannels, nChannels, p) - recurrence matrix (nChannels - number of signals,
%                                      p - model order)
%         nf - frequency resolution
%
%     Output:
%         AL(nf, nChannels, nChannels)
%     '''
new=1;
if new
    [nChannels, nChannels, p] = size(A);
    %# expoentes sao o array de expoentes da fft, com todas frequencias por todos lags
    %   exponents = reshape((-Jimag*pi*kron(0:(nf-1),(1:p))/nf),nf,p);
    exponents = reshape((-1i*pi*kron(0:(nf-1),(1:p))/nf),p,nf).';
    %    # Af faz a multiplicacoes da exp(ar) pela matrix A, para todas frequencias
    %    # as funcoes repeat e transpose sao truques para possibilitar o calculo vetorial
    Areshaped=reshape(A, nChannels,nChannels,1,p);
    
    %Af = (Areshaped.repeat(nf, axis=3))
    Af=zeros(nChannels,nChannels,nf,p);
    for kk=1:nf
        Af(:,:,kk,:)=Areshaped;
    end;
    for i=1:nChannels,
        for k=1:nChannels,
            Af(i,k,:,:)=reshape(Af(i,k,:,:),nf,p).*exp(exponents);
        end;
    end;
    
    Af=permute(Af, [3,1,2,4]);
    %    # o fft soma o valor para todos os lags
    AL=zeros(nChannels,nChannels,nf);
    H=zeros(nChannels,nChannels,nf);
    for kk=1:nf,
        temp=zeros(nChannels,nChannels);
        for k=1:p
            temp = temp+reshape(Af(kk,:,:,k),nChannels,nChannels);
        end;
        temp=eye(nChannels)-temp;
        AL(:,:,kk) = reshape(temp,1,nChannels,nChannels);
        H(:,:,kk)=inv(AL(:,:,kk));
    end;
else
    [K1,K2,p] = size(A);
    z = 1i*2*pi/Fs;
    for n=1:nf,
        atmp = zeros(K1);
        for k = 1:p+1,
            atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;
        AL(:,:,n)=atmp;
        H(:,:,n)  = inv(atmp);

    end
end