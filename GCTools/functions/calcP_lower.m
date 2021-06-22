function [P_lower]=calcP_lower(V_xyz,dim_xyz)

Sgma_xx = V_xyz(1:dim_xyz(1),1:dim_xyz(1));
Sgma_xy = V_xyz(1:dim_xyz(1),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2));
Sgma_xz = V_xyz(1:dim_xyz(1),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3));
Sgma_yx = V_xyz(dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),1:dim_xyz(1));
Sgma_yy = V_xyz(dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2));
Sgma_yz = V_xyz(dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3));
Sgma_zx = V_xyz(dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),1:dim_xyz(1));
Sgma_zy = V_xyz(dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),dim_xyz(1)+1:dim_xyz(1)+dim_xyz(2));
Sgma_zz = V_xyz(dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3),dim_xyz(1)+dim_xyz(2)+1:dim_xyz(1)+dim_xyz(2)+dim_xyz(3));

P1 = [eye(dim_xyz(1))     zeros(dim_xyz(1),dim_xyz(2))     zeros(dim_xyz(1),dim_xyz(3));
    -Sgma_yx/(Sgma_xx)    eye(dim_xyz(2))                  zeros(dim_xyz(2),dim_xyz(3));
    -Sgma_zx/(Sgma_xx)    zeros(dim_xyz(3),dim_xyz(2))     eye(dim_xyz(3))];

P2_temp = -(Sgma_zy - Sgma_zx/(Sgma_xx)*Sgma_xy)/(Sgma_yy - Sgma_yx/(Sgma_xx)*Sgma_xy);

P2 = [eye(dim_xyz(1))               zeros(dim_xyz(1),dim_xyz(2))         zeros(dim_xyz(1),dim_xyz(3));
    zeros(dim_xyz(2),dim_xyz(1))      eye(dim_xyz(2))                  zeros(dim_xyz(2),dim_xyz(3));
    zeros(dim_xyz(3),dim_xyz(1))           P2_temp                            eye(dim_xyz(3))];

P_lower = P2*P1;

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