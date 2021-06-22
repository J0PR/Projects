function tGC=newGCTime(X,a,covM)
% X (nvars*N)
% a (nvars*(nvars*m))
% covM (nvars*nvars)

% This method follows Eq. 19 from [1]
% Refs:
% [1] S. Hu et al., "Causality analysis of neural connectivity: critical examination of
% existing methods and advances of new methods.", IEEE transactions on neural networks, 2011

[nvars N]=size(X);
[nvars nvars_m]=size(a);
m=nvars_m/nvars;
%reshape (nvars*(nvars*P))->((nvars*nvars)*P)
a=reshape(a,[nvars nvars m]);

tGC=zeros(nvars,nvars);
for i=1:nvars
    for k=1:nvars
        if i~=k
            numerator = 0;
            for t=m+1:N
                temp_numerator=0;
                for j=1:m
                    temp_numerator = temp_numerator + a(k,i,j)*X(i,t-j);
                end
                numerator=numerator+temp_numerator^2;
            end
            
            den1=0;
            den2=0;
            for h=1:nvars
                if h~=k
                    for t=m+1:N
                        temp_den1=0;
                        for j=1:m
                            temp_den1 = temp_den1 + a(k,h,j)*X(h,t-j);
                        end
                        den1 = den1 + temp_den1^2;
                    end
                end
            end
            den2= N*covM(k,k);
%             residual array was replaced by covariance matrix.
%             for t=m+1:N-(m+1)
%                 den2 = den2 + (er(k,t))^2;
%             end
            tGC(k,i)=numerator/(den1+den2);
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