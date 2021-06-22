function [B,D]=extend_AR(A,covM)

[nvars nvars_p]=size(A);
p=nvars_p/nvars;

[L,D]=ldl(covM);

A=reshape(A,[nvars nvars p]);
B=zeros(size(A));
for i =1:p
    B(:,:,i)=L\A(:,:,i); %inv(L)*A(:,:,i);
end
%B0=eye(nvars)-inv(L); B0 is not used. Metrics use only causal part for
%now.
B=reshape(B,[nvars nvars*p]);