function [Sx, St, freqs] = xST(data, Fs)
%cross Stockwell Transform

[N M] = size(data);
for i=1:N
    for j=1:i
         [Si]=stran(data(i,:));
         [Sj]=stran(data(j,:));
        Sx(j,i,:,:)=Sj.*conj(Si);
        Sx(i,j,:,:)=Si.*conj(Sj);
    end
    St(i,:,:)=Si;%Stockwell tranform
end
freqs = linspace(0,Fs/2,size(St,2));