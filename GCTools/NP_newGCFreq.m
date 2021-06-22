function [fGC freqs] = NP_newGCFreq(X,Fs,NLAGS)
%UNDER CONSTRUCTION

[H Z S freqs]=calcHandCovNParam(X,Fs);

outDTF=DTF('np',H,freqs,Fs);
AR=leftDivision(eye(K1),H);
%converter para AR temporal e calcular o er.

fGC=newGCFreq(X,AR,er,Fs,freqs);


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