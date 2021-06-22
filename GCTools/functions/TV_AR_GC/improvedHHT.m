function [f, phi]=improvedHHT(y,Fs,last_phi)
%To use only in one IMF! 
% phase information in PHI is not reliable...
%y : nPoints*1
if size(y,2)>size(y,1)
    y=y';
end
initial_length=0.1*length(y); % % of data to include in initial data.
initial_data=y(1:initial_length,:);
strv=VARMA_est(initial_data,2,1);
[n, n, p]=size(strv.lastPhi(:,:,2:end));
ret=tvKalman(strv.lastPhi(:,:,2:end),strv.lastTheta(:,:,2:end),y');
lambdas=zeros(n*p,size(ret.A,3));
%f=zeros(n,size(ret.A,3));
for i=1:size(ret.A,3)
    lambdas(:,i)=eig(ret.A(:,:,i));
end

f=abs(log(lambdas(1,:)))/(2*pi*(1/Fs));
f(f>Fs/2)=Fs/2;
if nargin==3%calcular também as fases
    %começar a criar as fases atraves da ultima fase e ir subraidno as
    %frequencias
    f_rad=abs(log(lambdas(1,:)));
    f_diffs=cumsum(f_rad(end:-1:1));
    f_diffs=f_diffs(end:-1:1);
    phi=last_phi-f_diffs;
else
    phi=[];
end


