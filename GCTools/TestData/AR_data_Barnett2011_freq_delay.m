function [X,f,w]=AR_data_Barnett2011_freq_delay(N,Fs,f,delay,caus)
% Minimal VAR(1) from [1] (true model order 1)...
% According to [1] a=0.3, b=-0.8 and c=0.2104 for significant causality
% [1] Behaviour of Granger causality under filtering: theoretical invariance and practical application.
%a=0.3; b=-0.8;d=-0.99;%w=acos(b*(d-1)/(4*d));f=w*Fs/(2*pi);
a=0.3;  b=1.337025;d=-0.98;w=acos(b*(d-1)/(4*d));f=w*Fs/(2*pi);
%w=f*(2*pi)/Fs;

delay=delay*Fs;

if nargin ==1
     c=0.2104;
else
    % by solving this:  syms b caus c; solve(caus==log((1+b^2+c^2+sqrt((1+b^2+c^2)^2-4*b^2))/2), c);
    %updated for 'd': fYX_lambda = log(c^2/(b^2 + 2*cos(lambda)*b*d - 2*cos(lambda)*b + d^2 - 2*cos(2*lambda)*d + 1) + 1)
    % yields:
    lambda=w;
    c=((exp(caus) - 1)/(2*d - 2*b*cos(lambda) + b^2 + d^2 - 4*d*cos(lambda)^2 + 2*b*d*cos(lambda) + 1))^(1/2)*(b^2 + 2*cos(lambda)*b*d - 2*cos(lambda)*b + d^2 - 2*cos(2*lambda)*d + 1);
end

nstd1=1;
nstd2=1;

N=N+40+delay;
x=zeros(1,N+2);
y=zeros(1,N+2);

x(1) = 1;
y(1) = 1;
change=1;
for n=3+delay:N+2
    %y->x
    %     if n>=(N+40+1)/2 && change
    %         a=0.9;
    %         change=0;
    %         nstd1=1;
    %     end
    if caus == 0
        x(n)=b*x(n-1) + d*x(n-2) + c*y(n-delay) + normrnd(0,nstd1); %stdvar 1
    else
        x(n)=a*x(n-1) + c*y(n-delay) + normrnd(0,nstd1); %stdvar 1
    end
    y(n)=b*y(n-1) + d*y(n-2) + normrnd(0,nstd2); %stdvar 1

end


x=x(3+40+delay:N+2);
y=y(3+40+delay:N+2);

X=[y;x];
X=X./repmat(std(X,0,2),[1 size(X,2)]);
X=X-repmat(mean(X,2),[1 size(X,2)]);

%plot(X');
