function   X=glags(x,lags)
%
% this function generates a matrix containing lags of x
%
% Input arguments:
% x: an (n) x (1) vector series
% lags: the number of lags to be generated
%
% Output argument:
% X: an (n-lags) x (lags) matrix containing the lagged variables
% Note that lags observations are lost
%

[n,junk]=size(x);
X=[];
for i=1:lags
 X=[X x(lags-i+1:n-i,:)];
end
