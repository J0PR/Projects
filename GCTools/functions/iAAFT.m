function s=iAAFT(y,c)
% This Matlab version was adapted from Victor Venema's surrogates tools,
% Victor.Venema@uni-bonn.de, http:\\www.meteo.uni-bonn.de\victor, or 
% http:\\www.meteo.uni-bonn.de\victor\themes\surrogates\
% for the generation of surrogate cloud fields. 
% First version: May 2003.
% This version:  November 2003.

for i=1:c
    meanValue = mean(y);
    sorted_values = sort(y - meanValue);
    fourier_coeff = abs(ifft(y - meanValue))';
    [surrogate, errorAmplitude, errorSpec] = iaaft_loop_1d(fourier_coeff', sorted_values');
    surrogate = surrogate + meanValue;
    s(:,i)=surrogate;
end
        
function [y, errorAmplitude, errorSpec] = iaaft_loop_1d(fourierCoeff, sortedValues, initialisationStr, makePlots, template)
% INPUT: 
% fourierCoeff:   The 1 dimensional Fourier coefficients that describe the
%                 structure and implicitely pass the size of the matrix
% sortedValues:   A vector with all the wanted amplitudes (e.g. LWC of LWP
%                 values) sorted in acending order.
% OUTPUT:
% y:              The 1D IAAFT surrogate time series
% errorAmplitude: The amount of addaption that was made in the last
%                 amplitude addaption relative to the total standard deviation.
% errorSpec:      The amont of addaption that was made in the last fourier
%                 coefficient addaption relative to the total standard deviation

% Settings
errorThresshold = 2e-4; %
timeThresshold  = Inf;  % Time in seconds or Inf to remove this condition
speedThresshold = 1e-5; % Minimal convergence speed in the maximum error.
verbose = 1;
makePlots = 0; % Best used together with debugging.

% Initialse function
noValues = size(fourierCoeff);
errorAmplitude = 1;
errorSpec = 1;
oldTotalError = 100;
speed = 1;
standardDeviation = std(sortedValues);
t = cputime;

% The method starts with a randomized uncorrelated time series y with the pdf of
% sorted_values
[dummy,index]=sort(rand(size(sortedValues)));
y(index) = sortedValues;

% Main intative loop
while ( (errorAmplitude > errorThresshold | errorSpec > errorThresshold) & (cputime-t < timeThresshold) & (speed > speedThresshold) )
    % addapt the power spectrum
    oldSurrogate = y;    
    x=ifft(y);
    phase = angle(x);
    x = fourierCoeff .* exp(i*phase);
    y = fft(x);
    difference=mean(mean(abs(real(y)-real(oldSurrogate))));
    errorSpec = difference/standardDeviation;
    
    
    if (makePlots)
        plot(real(y))
        title('Surrogate after spectal adaptation')
        axis tight
        pause(0.01)
    end
        
    % addept the amplitude distribution
    oldSurrogate = y;
    [dummy,index]=sort(real(y));
    y(index)=sortedValues;
    difference=mean(mean(abs(real(y)-real(oldSurrogate))));
    errorAmplitude = difference/standardDeviation;
    
    
    if (makePlots)
        plot(real(y))
        title('Surrogate after amplitude adaptation')
        axis tight
        pause(0.01)
    end
    totalError = errorSpec + errorAmplitude;
    speed = (oldTotalError - totalError) / totalError;
    
    oldTotalError = totalError;
end

y = real(y);
