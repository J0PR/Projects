function [psi,freqs]=calcPSI(data,segleng,epleng,freqs,Fs,thres)

[Nsigs npoints]=size(data);
if isempty(segleng)
    segleng=200;
end
if segleng == 0;
    segleng = [];
end
if epleng == 0
    epleng= [];
end
% if isempty(epleng)
% epleng=floor(npoints/6);
% end
if segleng>npoints
    segleng=npoints;
end
nyq_bin=segleng/2+1;
nyq_freq=Fs/2;
fres=0.1;
if strcmp(freqs,'all')
    b1 = 1:nyq_bin-1;
    b2 = 2:nyq_bin;
    freqs=zeros(length(b1),2);
    freqs(:,1)=b1;
    freqs(:,2)=b2;
elseif strcmp(freqs,'mean')
    freqs = [];
else
    freqs=round(freqs*nyq_bin/nyq_freq);
    freqs(freqs==0)=1;
end

[psi, stdpsi, psisum, stdpsisum]=data2psi(data',segleng,epleng,freqs);

if ~isempty(epleng)
    psiVal = psi./(stdpsi+eps);
    if ~isempty(thres)
        psi(psiVal<thres)=0;
    else
        psi(psiVal<0)=0;
    end
else
    
    psi(psi<0)=0;

end
psi=permute(psi, [2 1 3]);

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



