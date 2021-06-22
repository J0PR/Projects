function [fGC freqs] = CP_newGCFreq(X,Fs,NLAGS)
freqs   =   [0:Fs/200:Fs/2-Fs/200];
cp_mvar = mvar_0lag_out(X,NLAGS,'CP');
AR=cp_mvar.alpha;
er=cp_mvar.u';
fGC=newGCFreq(X,AR,er,Fs,freqs);

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