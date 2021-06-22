function [H_temp Sgma_bar_f] = partitionImprov(H_f, V_xyz_norm,dim_xyz)
%Matrix Partition Improvement
%Calculate H_bar and Sgma_bar
H_temp = zeros(dim_xyz(1)+dim_xyz(3),dim_xyz(1)+dim_xyz(3),size(H_f,3));
H_bar = zeros(dim_xyz(1)+dim_xyz(3),dim_xyz(2),size(H_f,3));
Sgma_bar_f = zeros(dim_xyz(1)+dim_xyz(3),dim_xyz(1)+dim_xyz(3),size(H_f,3));
for ii = 1:size(H_f,3)
    H_f_xx = H_f(1:dim_xyz(1),1:dim_xyz(1),ii);
    H_f_xz = H_f(1:dim_xyz(1),dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),ii);
    H_f_zx = H_f(dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),1:dim_xyz(1),ii);
    H_f_zz = H_f(dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),ii);
    H_f_xy = H_f(1:dim_xyz(1),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),ii);
    H_f_zy = H_f(dim_xyz(1)+dim_xyz(2)+1:sum(dim_xyz),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),ii);
    
    H_temp(:,:,ii) = [H_f_xx H_f_xz;H_f_zx H_f_zz];
    H_bar(:,:,ii) = (H_temp(:,:,ii))\[H_f_xy;H_f_zy];   %for calculating Sgma_bar_f (2x2)
    Sgma_bar_f(:,:,ii) = [V_xyz_norm(1:dim_xyz(1),1:dim_xyz(1)) V_xyz_norm(1:dim_xyz(1),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3));...
        V_xyz_norm(dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),1:dim_xyz(1)) V_xyz_norm(dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3))] + ...
        H_bar(:,:,ii)*V_xyz_norm(dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2))*H_bar(:,:,ii)';
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