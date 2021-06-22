function [outgPDC,freqStruct] = LS_gPDC(X,Fs,NLAGS)

% frequency range to analyze (spectral analysis only)
freqs   =   [0:Fs/200:Fs/2-Fs/200];
% Estimate AR parameters
cp_mvar = mvar_0lag_out(X,NLAGS,'LS');
AR=cp_mvar.alpha;
PE=cp_mvar.Z;
M = size(AR,1);
% The PDF and the DTF can be displayed with the following functions
Xs.A = [eye(M),-AR]; Xs.B = eye(M); Xs.C  = PE; %Xs.C  = PE(:,(1-M:0)+end);
Xs.SampleRate = Fs;

outgPDC=gPDC(Xs.B,Xs.A,Xs.C,freqs,Fs);
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