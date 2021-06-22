function [GCsig IGC] = CP_GCTime(X,NLAGS,PVAL)

cpgc=tGC_0lag_out(X,NLAGS);

[PR,q] = bonferroni_significance(cpgc,PVAL);

GC = cpgc.gc;
GCsig = GC.*PR;
IGC=cpgc.igc;


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