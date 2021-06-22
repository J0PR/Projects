function [P_l P_u]=calculatePMatrix(C)

[N,N,F]=size(C);
C_l=C;
C_u=C;
if F==1
P_l=eye(N,N);
P_u=eye(N,N);
else
    P_l=repmat(eye(N,N),[1 1 F]);
    P_u=repmat(eye(N,N),[1 1 F]);
end
for k=1:F
    for j=1:N-1
        Pj_l=eye(N,N);
        for i=j+1:N
            Pj_l(i,j)=-C_l(i,j,k)/C_l(j,j,k);
        end
        C_l(:,:,k)=Pj_l*C_l(:,:,k);
        P_l(:,:,k)=Pj_l*P_l(:,:,k);
    end
end
for k=1:F
    for j=N:-1:2
        Pj_u=eye(N,N);
        for i=j-1:-1:1
            Pj_u(i,j)=-C_u(i,j,k)/C_u(j,j,k);
        end
        C_u(:,:,k)=Pj_u*C_u(:,:,k);
        P_u(:,:,k)=Pj_u*P_u(:,:,k);
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