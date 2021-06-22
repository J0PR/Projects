function [res auxStruct]=TV_ProcessData(data,optStruct,auxStruct)


measure=optStruct.measure;
type=optStruct.type;
MINP=optStruct.MINP;
MAXP=optStruct.MAXP;
info_crit=optStruct.info_crit;
pval=optStruct.pval;
nullPopMethod=optStruct.nullPopMethod;
numSurrogates=optStruct.numSurrogates;

[Nsigs npoints]=size(data);

%% detrend and demean
data = detrend(data')';
data=data-repmat(mean(data,2),[1 size(data,2)]);

%% Null population causality
if isempty(auxStruct)
    
    global glob_morder;
    if strcmp(type,'cond')
        if MINP == MAXP
            auxStruct.p=MAXP;
        else
            if isempty(glob_morder)
                auxStruct.p=model_order_cond(data,MINP,MAXP,optStruct);
            elseif glob_morder== -1
                auxStruct.p=model_order_cond(data,MINP,MAXP,optStruct);
                glob_morder=auxStruct.p;
            else
                auxStruct.p=glob_morder;
            end
        end
    else
        if MINP == MAXP
            auxStruct.p=MAXP*ones(Nsigs);
        else
            if isempty(glob_morder)
                auxStruct.p=model_order_pwise(data,MINP,MAXP,optStruct);
            elseif glob_morder== -1
                auxStruct.p=model_order_pwise(data,MINP,MAXP,optStruct);
                glob_morder=auxStruct.p;
            else
                auxStruct.p=glob_morder;
            end
        end
    end
    
    auxStruct.nullPop=buildNullPopulation(data,numSurrogates,nullPopMethod);
    nullRes=[];
    for i=1:size(auxStruct.nullPop,2)
    multiWaitbar( 'Calculating Null Dist', 0, 'Color', [0.2 0.6 0.2] );
        if strcmp(type,'cond')
            [nullRes{i}]=conditionalGCs(squeeze(auxStruct.nullPop(:,i,:)),optStruct,auxStruct.p);
            if isempty(nullRes{i})
                disp('Causality measure was not computed.')
                res=[];
                return;
            end
        else
            [nullRes{i}]=pairwiseGCs(squeeze(auxStruct.nullPop(:,i,:)),optStruct,auxStruct.p);
            if isempty(nullRes{i})
                disp('Causality measure was not computed.')
                res=[];
                return;
            end
        end
        multiWaitbar( 'Calculating Null Dist', i/size(auxStruct.nullPop,2));
    end
    multiWaitbar( 'Calculating Null Dist', 'Close');
    
    if strcmp(type,'cond')
        [res]=conditionalGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            disp('Causality measure was not computed.')
            return;
        end
    else
        [res]=pairwiseGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            disp('Causality measure was not computed.')
            return;
        end
    end
    res=TV_getSignificant(res,nullRes,pval,measure);
    
else
    %% sample causality
    if strcmp(type,'cond')
        [res]=conditionalGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            disp('Causality measure was not computed.')
            return;
        end
    else
        [res]=pairwiseGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            disp('Causality measure was not computed.')
            return;
        end
    end
end

%% DOI
if optStruct.DOI
    res.(measure)=DOI(res.(measure)); % calculate DOI
end

end

%% conditional casusality
function [res]=conditionalGCs(data,optStruct,p)
%data Nsigs*npoints
optStruct.NLAGS=p;
[res.(optStruct.measure), res.freqStruct] = TV_GrangerMetrics(data,optStruct);

end

%% pairwise causality
function res = pairwiseGCs(data,optStruct,pMatrix)
%data Nsigs*npoints

[Nsigs, npoints]=size(data);

pair_data=zeros(2,npoints);

for i=1:Nsigs
    for j=1:i
        if j~=i
            pair_data(1,:)=data(i,:);
            pair_data(2,:)=data(j,:);
            optStruct.NLAGS=pMatrix(i,j);
            [tempRes,freqStruct]=TV_GrangerMetrics(pair_data,optStruct);
            tempResTotal(i,j,:,:)=tempRes(1,2,:,:);
            tempResTotal(j,i,:,:)=tempRes(2,1,:,:);
        end
    end
end
res.(optStruct.measure)=tempResTotal;
res.freqStruct=freqStruct;
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

function ret=get_freqs(A,freqStruct,freqs) %deprecated.. 
%A=Nvars*Nvars*freqs*time
%freqs=Nbands*2*time
if strcmp(freqs,'all')
    ret=A;
elseif strcmp(freqs,'mean')
    ret=squeeze(mean(A,3));
elseif ndims(freqs)==3 %TV-freq_bands
    ret=zeros([size(A,1) size(A,2) size(freqs,1) size(A,4)]);
    for j=1:size(A,4)
        for i=1:size(freqs,1)
            ret(:,:,i,j)=squeeze(mean(A(:,:,freqStruct>=freqs(i,1,j) & freqStruct<freqs(i,2,j),j),3));
        end
    end
elseif ndims(freqs)==2 %TV-IFs for each driving var (each var is an IMF)
    ret=zeros([size(A,1) size(A,2) 1 size(A,4)]);
    for j=1:size(A,4)
        for k=1:size(freqs,1)
            ret(:,k,1,j)=squeeze(mean(A(:,k,freqStruct>=(freqs(k,j)-5) & freqStruct<(freqs(k,j)+5),j),3));
        end
    end
end
end
