function [bic,aic] = find_model_order(X,MINP,MAXP)
%adapted from CCA_Toolbox
[nvar,nobs] = size(X);

if(MAXP<=MINP) error('MAXP must be bigger than MINP, exiting'); end

bc = ones(1,MAXP).*999;
ac = ones(1,MAXP).*999;
for i = MINP:MAXP
    eval('res = cca_regress(X,i,0);','res = -1');   % estimate regression model, catch errors
    multiWaitbar( 'Finding model order', (i-MINP)/(MAXP-MINP));
    if(~isnumeric(res))
        [bc(i),ac(i)] = findbic(res,nvar,nobs,i);
        
    else
        disp('VAR failed');
        bc(i) = 999; 
        ac(i) = 999;
    end
end
multiWaitbar( 'Finding model order', 'Close');
[bicmin,bic] = min(bc);
[aicmin,aic] = min(ac);
end
%---------------------------------------------------------------------
function [bc,ac] = findbic(res,nvar,nobs,nlag)

error = log(det(res.Z));
nest = nvar*nvar*nlag;       
bc = error + (log(nobs)*nest/nobs);   
ac = error + (2*nest/nobs);
end


%    partially written by João Rodrigues
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