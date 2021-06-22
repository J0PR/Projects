function TEmat=gaussian_TE(data,d,tau,u)
%data Nsigs*npoints


[Nsigs npoints]=size(data);

d=3;
if d>=u;
    d=d-u+1;
end

TEmat = ones(Nsigs,Nsigs).*NaN;
for i=1:Nsigs
    for j=1:i
        if j~=i
            if Nsigs == 2
                conditional_data=[];
            else
                all_except_j_i = ~ismember(1:Nsigs, [j i]);
                conditional_data = data(all_except_j_i,:);
            end
            
            TEmat(i,j) = transferEntropy(data(j,:),data(i,:),conditional_data,d,tau,u);
            TEmat(j,i) = transferEntropy(data(i,:),data(j,:),conditional_data,d,tau,u);
        end
    end
end

end

function te=transferEntropy(y,x,z,d,tau,u)
%knn
for i=1:d
    xd(i,:)=x(1+d*tau-i*tau:end-i*tau-u);
end
if size(z,1) ~= 0
    for j=1:size(z,1)
        zd = z(j,:);
        for i=1:d
            xd(i+j*d,:)=zd(1+d*tau-(i-1)*tau:end-(i-1)*tau-u);
        end
    end
end
te=MI_gaussian(y(1+d*tau:end-u)',vertcat(x(1+u+d*tau:end),xd)')-MI_gaussian(y(1+d*tau:end-u)',xd');
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