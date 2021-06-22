function [res, auxStruct]=ProcessData(data,optStruct,auxStruct)

AR_mode=optStruct.AR_mode;
measure=optStruct.measure;
type=optStruct.type;
MINP=optStruct.MINP;
MAXP=optStruct.MAXP;
pval=optStruct.pval;
nullPopMethod=optStruct.nullPopMethod;
numSurrogates=optStruct.numSurrogates;

[Nsigs, npoints]=size(data);


%% detrend and demean
if strcmp(AR_mode,'NTrials_MVAR_Burg')
    data = cca_detrend_mtrial(data,optStruct.Nr,optStruct.Nl);
    data = cca_rm_temporalmean_mtrial(data,optStruct.Nr,optStruct.Nl,0);
elseif ~isfield(optStruct,'Nr') % in case there are repetitions and the AR_mode is not NTrials_MVAR_Burg, data must have been detrended previously
    data = detrend(data')';
    data=data-repmat(mean(data,2),[1 size(data,2)]);
end

%% PSI
if strcmp(optStruct.measure,'PSI')
    if ~isfield(optStruct,'segleng')
        optStruct.segleng=[];
    end
    if ~isfield(optStruct,'epleng')
        optStruct.epleng=[];
    end
    [res.(optStruct.measure),res.freqStruct] = calcPSI(data,optStruct.segleng,optStruct.epleng,optStruct.freqs,optStruct.Fs,optStruct.thres);
    return;
end

%% Null population causality
if isempty(auxStruct)
    
    if ~strcmp(measure,'TE')
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
    else
        auxStruct.p=[];
    end
    if ~strcmp(measure,'GCI')
        auxStruct.nullPop=buildNullPopulation(data,numSurrogates,nullPopMethod);
        nullRes=[];
        
        for i=1:size(auxStruct.nullPop,2)
            multiWaitbar( 'Calculating Null Dist', 0, 'Color', [0.2 0.6 0.2] );
            if strcmp(type,'cond')
                [nullRes{i}]=conditionalGCs(squeeze(auxStruct.nullPop(:,i,:)),optStruct,auxStruct.p);
                if isempty(nullRes{i})
                    error('Causality measure was not computed.')
                end
            else
                [nullRes{i}]=pairwiseGCs(squeeze(auxStruct.nullPop(:,i,:)),optStruct,auxStruct.p);
                if isempty(nullRes{i})
                    error('Causality measure was not computed.')
                end
            end
            multiWaitbar( 'Calculating Null Dist', i/size(auxStruct.nullPop,2));
        end
        multiWaitbar( 'Calculating Null Dist', 'Close');
    end
    if strcmp(type,'cond')
        [res]=conditionalGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            error('Causality measure was not computed.')
        end
    else
        [res]=pairwiseGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            error('Causality measure was not computed.')
        end
    end
    if ~strcmp(measure,'GCI')
        res=getSignificant(res,nullRes,pval,measure);
    end
else 
    %% sample causality
    if strcmp(type,'cond')
        [res]=conditionalGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            error('Causality measure was not computed.')
        end
    else
        [res]=pairwiseGCs(data,optStruct,auxStruct.p);
        if isempty(res)
            error('Causality measure was not computed.')
        end
    end
end

%% DOI
if optStruct.DOI
    res.(measure)=DOI(res.(measure)); % calculate DOI
end
%%
end


%% conditional casusality
function [res]=conditionalGCs(data,optStruct,p)
%data Nsigs*npoints

if strcmp(optStruct.measure, 'TE')
    switch optStruct.stats
        case {'knn'}
            res.(optStruct.measure) = knn_TE(data,optStruct.d,optStruct.tau,optStruct.u);
        case {'gaussian'}
            res.(optStruct.measure) = gaussian_TE(data,optStruct.d,optStruct.tau,optStruct.u);
        otherwise
            error(['Statistics with ' optStruct.stats ' does not exist for TE.']);
    end
else
    optStruct.NLAGS=p;
    if optStruct.extended_AR && ~optStruct.pcgc
        [M_lower,res.freqStruct]=GrangerMetrics(data,optStruct);
        [M_upper,res.freqStruct]=GrangerMetrics(data(end:-1:1,:),optStruct);
        temp_GC=zeros(size(M_upper));
        for i=1:size(M_upper,3)
            M_upper(:,:,i)=M_upper(:,:,i)';
            temp_GC(:,:,i)=temp_GC(:,:,i)+triu(M_upper(:,:,i),1);
            temp_GC(:,:,i)=temp_GC(:,:,i)+tril(M_lower(:,:,i),-1);
        end
        res.(optStruct.measure)=temp_GC;
    else
        [res.(optStruct.measure),res.freqStruct] = GrangerMetrics(data,optStruct);
    end
end

end
%% pairwise causality
function resOut = pairwiseGCs(data,optStruct,pMatrix)
%data Nsigs*npoints

[Nsigs npoints]=size(data);

freqStruct=[];
pair_data=zeros(2,npoints);

for i=1:Nsigs
    for j=1:i
        if j~=i
            pair_data(1,:)=data(i,:);
            pair_data(2,:)=data(j,:);
            optStruct.NLAGS=pMatrix(i,j);
            if strcmp(optStruct.measure, 'TE')
                switch optStruct.stats
                    case {'knn'}
                        tempRes = knn_TE(pair_data,optStruct.d,optStruct.tau,optStruct.u);
                    case {'gaussian'}
                        tempRes = gaussian_TE(pair_data,optStruct.d,optStruct.tau,optStruct.u);
                    otherwise
                        error(['Statistics with ' optStruct.stats ' does not exist for TE.']);
                end
            else
                [tempRes,freqStruct] = GrangerMetrics(pair_data,optStruct);
            end
            tempResTotal(i,j,:)=tempRes(1,2,:);
            tempResTotal(j,i,:)=tempRes(2,1,:);
        end
    end
end
res=mainDiagNaN(tempResTotal);
resOut.(optStruct.measure)=res;
if ~isempty(freqStruct)
    resOut.freqStruct=freqStruct;
end
end


%% aux functions
function M=mainDiagNaN(M)
N1=size(M,1);N2=size(M,2);N3=size(M,3);
aux=repmat(eye(N1),[1 1 N3]);
aux(find(aux))=NaN;
M=M+aux;
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