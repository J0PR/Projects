function [TGCList,optwindow]=Step_Wise_GC(dtimeSA,order, maximumTW, pval)
%add toolbox GUIwithavGC folder to path!!
[Nv,Ns]=size(dtimeSA);
Nr=1;
for j=1:Nv
    for k=1:Nv
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % estimating causality for each direction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j ~= k
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % optimally dividing time windows
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if k > j
                % temporal
                    % optimal time window for each pair of TS in the
                    % data set
                    clear optchangepoint
                    [optcoeff, opterror, optchangepoint, BIC, optexit] = ...
                        opttvCau(dtimeSA([j,k],:),Nr,Ns,order, 10, maximumTW);
                    %optchangepoint
                    for i = 1 : size(optchangepoint, 1)
                        optwindow{j,k}(i,:) = [optchangepoint(i,1),optchangepoint(i,2)];
                    end
                    if size(optchangepoint, 1) == 0
                        optwindow{j,k}(1,:) = [1,Ns];
                    end
            else
                optwindow{j,k} = optwindow{k,j};
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % causality estimation from kth region to jth region
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [TGCAvg(j,k), TpvalAvg(j,k), f, fc, TGCList{j,k}, timefc, TpvalList{j,k}] = ...
                TCau(dtimeSA(j,:),dtimeSA(k,:),Nr,Ns,order,optwindow{j,k}, 0);
        end
    end
end
TGCList=DOIandStatVal(TGCList,TpvalList,pval);

function TGCList=DOIandStatVal(TGCList,TpvalList,pval)

for i=1:size(TGCList,1)
    for j=1:i
        if j~=i
            tempGC1=TGCList{i,j};
            tempGC2=TGCList{j,i};
            aux1=tempGC1>tempGC2;
            aux2=tempGC2>tempGC1;
            tempGC1=tempGC1.*aux1;
            tempGC2=tempGC2.*aux2;
            
            tempVal1=TpvalList{i,j};
            tempVal2=TpvalList{j,i};
            aux1=tempVal1>pval;
            aux2=tempVal2>pval;
            tempGC1=tempGC1.*aux1;
            tempGC2=tempGC2.*aux2;
            TGCList{i,j}=tempGC1;
            TGCList{j,i}=tempGC2;
        end
    end
end
