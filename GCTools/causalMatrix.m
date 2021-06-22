function [mres,freqs,p,sig_thresh]=causalMatrix(data,optStruct)
%adjMatrix=causalMatrix(data(nVars,nCols),optStruct);
optStruct.pcgc=0;
[res,auxStruct]=ProcessData(data,optStruct,[]);
if ~isempty(res)
    mres=res.(optStruct.measure);
    if isfield(auxStruct,'p')
        p=auxStruct.p;
    else
        p=[];
    end
    if isfield(res,'freqStruct')
        freqs=res.freqStruct;
    else
        freqs=[];
    end
    if isfield(res,'sig_thresh')
        sig_thresh=res.sig_thresh;
    else
        sig_thresh=[];
    end
else
    mres=[];
    freqs=[];
    p=[];
    sig_thresh=[];
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