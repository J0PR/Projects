function [outPDC,freqStruct] = MTrial_PDC(X,Nr,Nl,Fs,NLAGS,freqs)

% frequency range to analyze (spectral analysis only)
if nargin <6
    freqs   =   [0:Fs/200:Fs/2-Fs/200];
end
% Estimate AR parameters
[AR,PE]=armorf(X,Nr,Nl,NLAGS);
M = size(AR,1);
% The PDF and the DTF can be displayed with the following functions
Xs.A = [eye(M),-AR]; Xs.B = eye(M); Xs.C  = PE;
Xs.datatype = 'MVAR';
Xs.SampleRate = Fs;

outPDC=PDC(Xs.B,Xs.A,freqs,Fs);
freqStruct=freqs;

%    Written by João Rodrigues
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.