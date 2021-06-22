function s=AAFT(x,c)
%Syntax: s=AAFT(x,c)
%___________________
%
% Makes c Amplitude Adjusted Fourier Transformed (AAFT) surrogates of a time
% series x.
%
% s is the AAFT time series.
% x is the original time series.
% c is the number of surrogates.
%

if nargin<1 | isempty(x)==1
   error('You should provide a time series.');
else
   % x must be a vector
   if min(size(x))>1
      error('Invalid time series.');
   end
   x=x(:);
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

for i=1:c
    % Initialize
    y=x;
    % Make n normal random devaiates
    normal=sort(randn(size(y)));
    % Sort y and extract the ranks
    [y,T]=sort(y);
    [T,r]=sort(T);
    % Assign the ranks of y to the normal deviates and apply the phase
    %  randomization
    normal=phaseran(normal(r));
    % Extract the ranks of the phase randomized normal deviates
    [normal,T]=sort(normal);
    [T,r]=sort(T);
    % Assign the ranks of the phase randomized normal deviates to y and
    %  obtain the AAFT surrogates
    s(:,i)=y(r);
end