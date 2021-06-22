function res=getSignificant(res,nullRes,pval,measure)
%significance testing
if isempty(nullRes)
    return;
end
nullMat=zeros([size(res.(measure)) length(nullRes)]);
for i=1:length(nullRes)
    auxRes=nullRes{i};
    if ndims(auxRes.(measure)) == 2
        aux=auxRes.(measure);
        aux(find(eye(size(aux,1))))=NaN;
        nullMat(:,:,i)=aux;
    else
        aux=auxRes.(measure);
        aux(find(repmat(eye(size(aux,1)),[1 1 size(aux,3)])))=NaN;
        nullMat(:,:,:,i)=aux;
    end
end
if ndims(nullMat) == 3
    nullVals=nullMat(:);
    clear nullMat;
    nullVals(isnan(nullVals))=[];
    thres=prctile(nullVals,100-100*pval);
else
    for i=1:size(nullMat,3)
        aux=nullMat(:,:,i,:);
        for ii=1:size(aux,1)
            for jj=1:size(aux,2)
                nullVals=squeeze(aux(ii,jj,:,:));
                nullVals(isnan(nullVals))=[];
                thres(ii,jj,i)=prctile(nullVals,100-100*pval);
            end
        end
        clear aux;
    end
end

mat=res.(measure);
if ndims(mat) == 2
    mat(mat<thres)=0;
else
    for i=1:size(mat,3)
        tempMat=mat(:,:,i);
        tempMat(tempMat<thres(:,:,i))=0;
        mat(:,:,i)=tempMat;
    end
end
res.(measure)=mat;
res.sig_thresh=thres;
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