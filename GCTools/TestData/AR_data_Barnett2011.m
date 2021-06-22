function [X]=AR_data_Barnett2011(N,caus)
% Minimal VAR(1) from [1] (true model order 1)...
% According to [1] a=0.3, b=-0.8 and c=0.2104 for significant causality
% [1] Behaviour of Granger causality under filtering: theoretical invariance and practical application.
a=0.3; b=-0.8;
if nargin ==1
     c=0.2104;
else
    % by solving this:  syms b caus c; solve(caus==log((1+b^2+c^2+sqrt((1+b^2+c^2)^2-4*b^2))/2), c);
    % yields:
    c=exp(-caus/2)*((exp(caus) - 1)*(- b^2 + exp(caus)))^(1/2);
end

nstd1=1;
nstd2=1;

N=N+40;
x=zeros(1,N+1);
y=zeros(1,N+1);

x(1) = 1;
y(1) = 1;
change=1;
for n=2:N+1
    %y->x
%     if n>=(N+40+1)/2 && change
%         a=0.9;
%         change=0;
%         nstd1=1;
%     end
    x(n)=a*x(n-1) + c*y(n-1) + normrnd(0,nstd1); %stdvar 1
    y(n)=b*y(n-1) + normrnd(0,nstd2); %stdvar 1
end


x=x(2+40:N+1);
y=y(2+40:N+1);

X=[x;y];

%plot(X');
