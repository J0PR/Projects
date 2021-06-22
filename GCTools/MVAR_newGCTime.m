function tGC = MVAR_newGCTime(X,NLAGS)
[nvar,nobs] = size(X);
m = mean(X');
if(abs(sum(m)) > 0.0001)
    mall = repmat(m',1,nobs);
    X = X-mall;
end

[AR,RCF,PE,DC,u] = mvar(X', NLAGS);
er=u';
tGC=newGCTime(X,AR,er);

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