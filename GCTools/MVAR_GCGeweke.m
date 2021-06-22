function [outGGC,freqStruct] = MVAR_GCGeweke(X,Fs,NLAGS)

% frequency range to analyze (spectral analysis only)
freqs   =   [0:Fs/200:Fs/2-Fs/200];
% Estimate AR parameters
[AR,RC,PE] = mvar(X',NLAGS);
M = size(AR,1);
% The PDF and the DTF can be displayed with the following functions
Xs.A = [eye(M),-AR]; Xs.B = eye(M); Xs.C  = PE(:,(1-M:0)+end);
Xs.datatype = 'MVAR';
Xs.SampleRate = Fs;

outGGC=GewekeGC(X,Xs.B,Xs.A,Xs.C,NLAGS,freqs,Fs);
freqStruct=freqs;


%    Written by Jo�o Rodrigues
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