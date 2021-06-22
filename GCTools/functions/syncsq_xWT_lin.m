function [Wx, W, T, w, freqs] = syncsq_xWT_lin(data, Fs)

t=[0:size(data,2)-1]*(1/Fs);
nv=32;
[N M] = size(data);
for i=1:N
    for j=1:i
         [Ti, freqs, Wxi, as, wi]=synsq_cwt_fw_lin(t,data(i,:),nv);
         [Tj, freqs, Wxj, as, wj]=synsq_cwt_fw_lin(t,data(j,:),nv);
        Wx(j,i,:,:)=Wxj.*conj(Wxi);
        Wx(i,j,:,:)=Wxi.*conj(Wxj);
    end
    w(i,:,:)=wi;%inst freq map
    W(i,:,:)=Wxi;%wavelet tranform
    T(i,:,:)=Ti;%synchrosqueezed wavelet transform
end