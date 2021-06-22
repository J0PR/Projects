function [GCout IGC] = LS_GCTime(X,NLAGS,PVAL)

ret = cca_granger_regress(X,NLAGS,1);

[PR,q] = cca_findsignificance(ret,PVAL,1);

GC = ret.gc;
if PVAL ~=0 && PVAL ~=1
    GCout = GC.*PR;
else
    GCout = GC;
end
IGC=ret.igc;

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