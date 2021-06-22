function pMatrix=model_order_pwise(data,MINP,MAXP,optStruct)


info_crit=optStruct.info_crit;
[Nsigs npoints]=size(data);


pair_data=zeros(2,npoints);

for i=1:Nsigs
    for j=1:i
        if j~=i
            pair_data(1,:)=data(i,:);
            pair_data(2,:)=data(j,:);
            if strcmp(optStruct.AR_mode,'NTrials_MVAR_Burg')
                [bic,aic] = cca_find_model_order_mtrial(pair_data,optStruct.Nr,optStruct.Nl,MINP,MAXP);
            else
                [bic,aic] = find_model_order(pair_data,MINP,MAXP);
            end
            if strcmp(info_crit,'bic')
                pMatrix(i,j)=bic;
            else
                pMatrix(i,j)=aic;
            end
        end
    end
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