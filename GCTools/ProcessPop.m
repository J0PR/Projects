function [measuresStruct]=ProcessPop(sigsArray,measure,ARmode,processType,Fs,freqs,MINP,MAXP,pvalue,bypassSurr)
%freqs: varias bandas de frequencia como na psi. cada linha é uma banda.
multiWaitbar( 'Processing Population', 0, 'Color', [0.2 0.6 0.2] );

measuresStruct=cell(1,size(sigsArray,2));

for i=1:size(sigsArray,2)
    if i == 1 || bypassSurr == 0
        [msr nullPopStruct]=ProcessData(sigsArray{i},measure,ARmode,processType,Fs,freqs,MINP,MAXP,pvalue,[]);
    else
        [msr nullPopStruct]=ProcessData(sigsArray{i},measure,ARmode,processType,Fs,freqs,MINP,MAXP,pvalue,nullPopStruct);
    end
    measuresStruct{i}=msr.(measure);
    multiWaitbar( 'Processing Population', i/size(sigsArray,2));
end
multiWaitbar( 'Processing Population', 'Close');
end

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