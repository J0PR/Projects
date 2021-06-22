function [y ind]=info_gain(drive,X,nvar,ndmax)
indici=setdiff(1:nvar,drive); %eliminate the candidate driver from the set
t=X{drive};
Zt=[];
iter=0;
reverseStr = '';
 fprintf(1,'info gain:  ');
 s = clock;
for nd=1:ndmax
    n1=length(indici);
    z=zeros(n1,1);
    for k=1:n1
        iter=iter+1;
        Zd=[Zt X{indici(k)}];
        z(k)= MI_gaussian(Zd,t); %compute conditional MI, here from covariance matrix, but you can use the exact formula
        %printing and estimated time left
        is = etime(clock,s);
        esttime = (is/iter) * (ndmax*n1);
        msg = sprintf('%3.2f ETL: %3.2f',iter / (ndmax*n1),(esttime-etime(clock,s))/60);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        %fprintf(1,'\b%3.2f ETL: %3.2f',iter / (ndmax*n1),esttime-etime(clock,s));
        %printing and estimated time left
    end
    [y(1,nd) id]=max(z); %greedy algorithm, find the max contribution, store it and remove it from the set of candidates
    Zt=[Zt X{indici(id)}];
    ind(1,nd)=indici(id);
    indici=setdiff(indici,indici(id));
end
fprintf('\n')

function ret=MI(Zd,t)


ret = MIknn(Zd',t');