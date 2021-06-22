function res=TV_getSignificant(res,nullRes,pval,measure)
%significance testing
if isempty(nullRes)
    return;
end
timeDim=ndims(res.(measure));%time dimension is always the last dim.
for t=1:size(res.(measure),timeDim)
    if timeDim==4
        tRes=res.(measure)(:,:,:,t);
    else
        tRes=res.(measure)(:,:,t);
    end
    nullMat=zeros([size(tRes) length(nullRes)]);
    for i=1:length(nullRes)
        auxRes=nullRes{i};
        if timeDim==4
            tNullRes=auxRes.(measure)(:,:,:,t);
            tNullRes(find(repmat(eye(size(tNullRes,1)),[1 1 size(tNullRes,1)])))=NaN;
            nullMat(:,:,:,i)=tNullRes;
        else
            tNullRes=auxRes.(measure)(:,:,t);
            tNullRes(find(eye(size(tNullRes,1))))=NaN;
            nullMat(:,:,i)=tNullRes;
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
            nullVals=aux(:);
            clear aux;
            nullVals(isnan(nullVals))=[];
            thres(i)=prctile(nullVals,100-100*pval);
        end
    end
    
    mat=res.(measure);
    if ndims(mat) == 2
        mat(mat<thres)=0;
    else
        for i=1:size(mat,3)
            tempMat=mat(:,:,i);
            tempMat(tempMat<thres(i))=0;
            mat(:,:,i)=tempMat;
        end
    end
    if timeDim==4
        res.(measure)(:,:,:,t)=mat;
    else
        res.(measure)(:,:,t)=mat;
    end
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