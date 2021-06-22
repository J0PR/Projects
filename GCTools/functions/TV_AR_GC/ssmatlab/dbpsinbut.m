function [compbp,ferror]=dbpsinbut(D,Omegap1,Omegap2,Omegas2,Di,Thetac,Lambda)
% *************************************************************************  
% This function obtains the band-pass filter based on the Butterworth 
% tangent filter corresponding to the parameters D(1), D(2), Omegap1, 
% Omegap2, Omegas2. See "The Use of Butterworth Filters for Trend and Cycle
% Estimation in Economic Time Series", Gómez, V. (2001), Journal of 
% Business and Economic Statistics, 19, 365-373.
% The filter model is
%
%           z_t = s_t + n_t,
%   Alpha(B)s_t = num(B)b_t,    Var(b_t)=1
%
% where Alpha(z) = (1 - 2*Alph*z + z^2)^Di, num(z) = (1 - Alph*z)^Di, and 
% n_t and b_t are independent white noises.
% 
% The filter numerator is num*(1/sa). The filter denominator is den. Thus,
% 
% H(z) = (1/sa)*(num(z)/den(z)) 
%
% The other filter is
%
% G(z) = (sqrt(Lambda)/sa)*(Alpha(z)/den(z))
%
% Input parameters:
%     D       : a (1 x 2) array containing the design tolerances D1 and D2.
%               It can be empty.
%     Omegap1 : a number, design frequency Omegap1 divided by pi. Required.
%     Omegap2 : a number, design frequency Omegap2 divided by pi. Required.
%     Omegas2 : a number, design frequency Omegas2 divided by pi. It can be 
%               empty 
%     Di      : a number, the exponent in Alpha(z) and num(z). It can be
%               empty.
%     Thetac  : a number, the frequency, divided by pi, of gain .5 in the 
%               But. tan. filter. It can be empty.
%     Lambda  : a number, the signal to noise ratio (sigma^2_n/sigma^2_b)
%               in the But. tangent filter. It can be empty.
% Note: The usual specification is D, Omegap1, Omegap2 and Omegas2 (Di, 
%       Thetac and Lambda empty). Alternatively, the user can enter 
%       Omegap1, Omegap2, Di and Thetac (D, Omegas2 and Lambda empty) or 
%       Omegap1, Omegap2, Di and Lambda (D, Omegas2 and Thetac empty).
%
% Output parameters: compbp, a structure containing the following fields
%     .num    : a polynomial of degree Di, (1 - Alph*z)^Di
%     .den    : a polynomial of degree 2*Di
%     .sa     : a positive number (sigmaa)
%     .Alpha  : a polynomial of degree 2*Di, (1 - 2*Alph*z + z^2)^Di
%     .Alph   : a number, in Alpha(z) = (1 - 2*Alph*z + z^2)^Di
%     .Di     : a positive integer
%     .Thetac : a number, the frequency of gain .5 in the But. tan. filter
%     .Lambda : a positive number, the square root of the noise to signal 
%               ratio (sigma^2_n/sigma^2_b) in the But. tangent filter.
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos, 
% Subdireccion Gral. de Analisis y P.E., 
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should 
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

ferror=0;
if isempty(Omegap1) || isempty(Omegap2)
 disp('Omegap1 and Omegap2 must be entered in dbpsinbut')
 ferror=1;
 return
end
Omegap1 = Omegap1*pi;
Omegap2 = Omegap2*pi;

if ( nargin == 6 )  
 %specification Omegap1, Omegap2, Di and Thetac
 if isempty(Di) || isempty(Thetac)
  disp('Di and Thetac must be entered in dbpsinbut when there are')
  disp('six arguments')
  ferror=2;
  return
 end
 Thetac=Thetac*pi;
 dd = double(Di);
 Lambda = 1/((2*sin(Thetac/2))^Di);
elseif ( nargin == 7 ) 
 %specification Omegap1, Omegap2, Di and Lambda
 if isempty(Di) || isempty(Lambda)
  disp('Di and Lambda must be entered in dbpsinbut when there are')
  disp('seven arguments')
  ferror=3;
  return
 end
 dd = double(Di);
 Thetac = 2*asin(1/(2*exp(log(Lambda)/(dd+dd))));
 Lambda = sqrt(Lambda);
else
 %specification D, Omegap1, Omegap2 and Omegas2 
 if isempty(D) || isempty(Omegas2)
  disp('D and Omegas2 must be entered in dbpsinbut when there are less')
  disp('than six arguments')
  return
 end
 Omegas2 = Omegas2*pi;
% transformation in the frequency domain: THETA = OMEGA-OMEGA_P1
 thetap = Omegap2 - Omegap1;
 thetas = Omegas2 - Omegap1;
 
 % Find Di
 sum1 = log(D(1)) - log(1-D(1));
 sum2 = log(D(2)) - log(1-D(2));
 sum = sum1 + sum2;
 
 den1 = ((tan(thetap/2))^2)/(1 + (tan(thetap/2))^2);
 lden1 = log(den1);
 den2 = ((tan(thetas/2))^2)/(1 + (tan(thetas/2))^2);
 lden2 = log(den2);
 deno = lden1 - lden2;

 dd = sum/deno;
 Di = round(dd);
 dd = double(Di);
 
 % Find Thetac given Di
 sum = exp(sum1/dd);
 sum = sum/den1;
 sum = sqrt(1/(sum - 1));
 Thetac = 2*atan(sum);
 
 Lambda = 1/((2*sin(Thetac/2))^Di);
end 

% now transformation in the time domain
Alph = cos((Omegap1+Omegap2)/2)/cos((Omegap2-Omegap1)/2);
nterm2 = floor(Di/2);
if ( mod(Di,2) == 0 ) 
 nterm1 = 0;
else
 nterm1 = 1;
end
sa = Lambda;

% set up moving average polynomial for the aggregate series
b = (sin(Thetac/2))^2;

% derive the denominator
den = 1;
for i = 1:nterm2
 a = 1 - 2*b*cos((pi+2*(i-1)*pi)/dd)+b^2;
 alp0 = b + sqrt(a) + sqrt((b+sqrt(a))^2 - 1);
 alp2 = b + sqrt(a) - sqrt((b+sqrt(a))^2 - 1); 
 alp1 = 2*(b - sqrt(a));   
    
 sa = sa*abs(alp0);
 
 delta(1) = 1;
 alp1 = alp1/alp0;
 alp2 = alp2/alp0;
 
% do transformation
 delta(2) = Alph*(alp1-2);
 delta(3) = Alph*Alph*(1-alp1+alp2) - alp1;
 delta(4) = Alph*(alp1-2*alp2);
 delta(5) = alp2;
 delta=fliplr(delta);
 den=conv(delta,den);
 clear delta;
end

if ( nterm1 == 1 )  
 alp0 = sqrt(b) + sqrt(b + 1);
 alp1 = sqrt(b) - sqrt(b + 1);
 sa = sa*abs(alp0);
 alp1 = alp1/alp0;
 den=conv([-alp1 Alph*(alp1-1) 1],den);
end

% derive the numerator
delta = [-Alph 1];
num(1) = 1;
for i = 1:Di
 num=conv(delta,num);
end

% derive Alpha
delta = [1 -2*Alph 1];
Alpha(1) = 1;
for i = 1:Di
 Alpha=conv(delta,Alpha);
end

compbp.num=num;          % filter numerator
compbp.den=den;          % filter denominator
compbp.sa=sa;            % filter sa 
compbp.Alpha=Alpha;      % filter Alpha
compbp.Di=Di;            % filter Di
compbp.Lambda=Lambda;    % filter Lambda
compbp.Thetac=Thetac;    % filter Thetac
compbp.Alph=Alph;        % filter Alph

end
