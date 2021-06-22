function improvedHHT_DEKF(y,Fs,p)
if size(y,1)>size(y,2)
    y=y';
end
n=size(y,1);
I=eye(n);
res = DEKF(y,p);
phi=res.phi;

for k=1:size(phi,3)
    auxA=eye(p-1);
    auxA=[auxA zeros(p-1,1)];
    auxA=kron(auxA,I);
    A(:,:,k)=[phi(:,:,k);auxA];
end
f_temp=zeros(n*p,size(A,3));
%f=zeros(n,size(ret.A,3));
for i=1:size(A,3)
    lambdas(:,i)=eig(A(:,:,i));
end
% test1=(real(lambdas(1:2:end-1,:)) == real(lambdas(2:2:end,:)));
% test2=abs(imag(lambdas(1:2:end-1,:))) == abs(imag(lambdas(2:2:end,:)));
% test3=+(test1&test2);%+converts logical to double...
lambdaSingle=lambdas(1:2:end-1,:);
% f=abs(log(lambdaSingle.*test3))/(2*pi*(1/Fs));
f=abs(log(lambdaSingle))/(2*pi*(1/Fs));
f(f==Inf)=0;
f(f>Fs/2)=0;
figure, plot(f')