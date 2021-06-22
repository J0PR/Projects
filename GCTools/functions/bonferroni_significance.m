function [PR,q] = bonferroni_significance(ret,pval)
%-----------------------------------------------------------------------
% INPUTS:   ret:    output of cca_granger_regress.m (or similar)
%           pval:   desired P-value threshold
%
% OUTPUT:   PR:    nvar by nvar matrix of 1s (sig connections)
%           q:     threshold used (=pval for no correction) 
%
%-----------------------------------------------------------------------

prb = ret.prb;
nvar = length(prb);
PR = zeros(nvar);

% Bonferroni correction
    
    q = pval./(nvar*(nvar-1));
    PR(prb<q) = 1;
    
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