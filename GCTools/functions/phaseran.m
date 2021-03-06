function s=phaseran(x,c)
%Syntax: s=phaseran(x,c)
%_______________________
%
% Makes c phase randomized surrogates of a time series x.
%
% s is the phase randomized time series.
% x is the original time series.
% c is the number of surrogates.


if nargin<1 | isempty(x)==1
   error('You should provide a time series.');
else
   % x must be a vector
   if min(size(x))>1
      error('Invalid time series.');
   end
   x=x(:);
   % N is the time series length
   N=length(x);
end

if nargin<2 | isempty(c)==1
   c=1;
else
   % c must be scalar
   if sum(size(c))>2
      error('c must be scalar.');
   end
   % c must be greater or equal than 1
   if c<1
      error('c must be greater or equal than 1.');
   end
end


% FFT on x
y=fft(x);
% Magnitudes
m=abs(y);
% Angles
p=angle(y);

% Half of the data points
h=floor(N/2);

for j=1:c
   % Randomized phases
   if rem(N,2)==0
      p1=rand(h-1,1)*2*pi;
      p(2:N)=[p1' p(h+1) -flipud(p1)'];
      % Adjust the magnitudes
      m=[m(1:h+1);flipud(m(2:h))];
   else
      p1=rand(h,1)*2*pi;
      p(2:N)=[p1 -flipud(p1)];
   end
   % Back to the complex numbers
   s(:,j)=m.*exp(i*p);
   % Back to the time series (phase randomized surrogates)
   s(:,j)=real(ifft(s(:,j)));
end