function [DTF,PDC,DC,GGC,GCGeweke,PartitionsGCGeweke,freqStruct] = CP_CausalMeasures(X,Fs,NLAGS)

% frequency range to analyze (spectral analysis only)
freqs   =   [0:Fs/200:Fs/2-Fs/200];
multiWaitbar( 'Calculating measures.', 0, 'Color', [0.2 0.6 0.2] );
% Estimate AR parameters
cp_mvar = mvar_0lag_out(X,NLAGS,'CP');
AR=cp_mvar.alpha;
PE=cp_mvar.Z;
M = size(AR,1);
% The PDF and the DTF can be displayed with the following functions
Xs.A = [eye(M),-AR]; Xs.B = eye(M); Xs.C  = PE; %Xs.C  = PE(:,(1-M:0)+end);
Xs.datatype = 'MVAR';
Xs.SampleRate = Fs;

%GCtest=freqGCtest(X,Xs.B,Xs.A,Xs.C,NLAGS,freqs,Fs);
multiWaitbar( 'Calculating measures.', 0.1);
GCGeweke=CPGewekeGC(X,Xs.B,Xs.A,Xs.C,NLAGS,freqs,Fs);
multiWaitbar( 'Calculating measures.', 0.2);
PartitionsGCGeweke=GewekeGC1Partitions(X,Xs.B,Xs.A,Xs.C,NLAGS,freqs,Fs);
multiWaitbar( 'Calculating measures.', 0.3);
%PartitionsGCPCA=PCAGCPartitions(X,Xs.B,Xs.A,Xs.C,NLAGS,freqs,Fs);
multiWaitbar( 'Calculating measures.', 0.4);
[S,h,PDC,DTF,DC,GGC]=mvfreqz(Xs.B,Xs.A,Xs.C,freqs,Fs);
multiWaitbar( 'Calculating measures.', 0.5);
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