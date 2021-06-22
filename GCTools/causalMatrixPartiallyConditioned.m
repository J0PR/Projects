function [mres,freqs,p]=causalMatrixPartiallyConditioned(data,optStruct)

measure=optStruct.measure;
freqs=optStruct.freqs;
Fs=optStruct.Fs;
MINP=optStruct.MINP;
MAXP=optStruct.MAXP;
optStruct.type='cond';
optStruct.pcgc=1;
ndmax=optStruct.ndmax;

morder=model_order_cond(data,MINP,MAXP,optStruct);
if ~strcmp(freqs,'all') && ~strcmp(freqs,'mean') %freq bands.. have to filter for init_partial_conditioning_par
    if size(freqs,1)>1
        disp('Sorry, partially conditioned analysis cannot be performed for multiple bands. Please call partiallyConditionedCausalMatrix once for each band.');
        res=[];
        return;
    end
    fdata = eegfilt(double(data),Fs,freqs(1),freqs(2));
    [y, ind]=init_partial_conditioning_par(fdata',ndmax,morder);
else
    [y, ind]=init_partial_conditioning_par(data',ndmax,morder);
end
figure,plot(1:ndmax-1,diff(y'));
xlabel('n_d')
ylabel('\Deltay')
title('choose the point in x corresponding to the "elbow" of the plots')
nd=input('Number of conditional variables according to the graph: ');

nvar = size(data,1);
for drive=1:nvar
    for target=1:nvar
        if drive ~= target
            A=ind(drive,:);
            conz = A(~ismembc(A(:), target));
            conz = conz(1:nd);
            temp_data=[data(drive,:);data(target,:);data(conz,:)];
            
            [tempRes,nullStruct] = ProcessData(temp_data,optStruct,[]);
            if ~isempty(tempRes)
                tempMres=tempRes.(measure);
                mres(target,drive,:)=tempMres(2,1,:);
                if isfield(nullStruct,'p')
                    p(target,drive)=nullStruct.p;
                else
                    p=[];
                end
                if isfield(tempRes,'freqStruct')
                    freqs=tempRes.freqStruct;
                else
                    freqs=[];
                end
            else
                mres=[];
                freqs=[];
                p=[];
                return
            end
        end
    end
end


