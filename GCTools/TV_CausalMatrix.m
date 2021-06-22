function [mres, freqs]=TV_CausalMatrix(data,optStruct)
%adjMatrix=TV_CausalMatrix(data(nVars,nCols),'gPDC','TV_MVAR','cond',250,[5 30],[],'ST',1,30,'bic',0.05,100,[]);
mres=[];
[res nullStruct]=TV_ProcessData(data,optStruct,[]);
if ~isempty(res)
    mres=res.(optStruct.measure);
    freqs=res.freqStruct;
end
multiWaitbar( 'Calculating measures.', 'Close');

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