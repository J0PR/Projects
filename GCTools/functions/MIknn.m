function mi=MIknn(x,y,kneig)

% Calculate MI value between 2 vector of any dimension
% x = input data (channel*points)
% y = input data (channel*points)
% kneig = k nearest neigbor for MI algorithm


%default-values
if exist('kneig','var') ~=1
    kneig=6; 
end


% check input data if format is correct
[Ndx,Nx]=size(x);
if Ndx>Nx
    x=x';
    [Ndx,Nx]=size(x);
end
[Ndy,Ny]=size(y);
if Ndy>Ny
    y=y';
    [Ndy,Ny]=size(y);
end

if Nx~=Ny
    if Nx>Ny
        N=Ny;
    else
        N=Nx;
    end
    fprintf('The two input vectors must have the same length');
    fprintf('Caluculation using the %d datapoints',N);
    
else
    N=Nx;
end

zwsp=[x;y];

fast_func=exist('MI','file');

if fast_func == 3 %.mex exists
    res=MI(zwsp,Ndx,Ndy,N,kneig);
else
    res=MI_knn(zwsp,Ndx,Ndy,N,kneig);
    %or else: mex -g MI.cpp
    %call MI now
end

mi=res;


function ret=MI_knn(x,dimx,dimy,N,K)

x=x+rand(size(x))*1e-8;
x = zscore(x')';
min_array=min(x,[],2)';
max_array=max(x,[],2)';
x=x-repmat(min(x,[],2),[1 size(x,2)]);

psi=zeros([1 N+1]);
psi(2)=-0.57721566490153;
for i=2:N
    psi(i+1)=psi(i)+1/(i-1);
end

BOX1=N-5;
scal=BOX1./(max_array-min_array);

ret=mir_xnyn(x,dimx,dimy,N,K,psi,scal);

function mir=mir_xnyn(x, dimx, dimy, N, K,psi,scal)

nn=ones([1 K+1]);

phi(2:K+1)=psi(2:K+1)-1./(1:K);


BOX=1;
while 0.5*BOX*BOX*K<N
    BOX=BOX*2;
end
epsilon=4.0/BOX;
inveps=BOX/4;
BOX1=N-5;

if dimx>1
    xx(1:dimx,1:N)=x(1:dimx,1:N);
end
if dimy>1
    yy(1:dimy,1:N)=x(dimx+1:end,1:N);
end

[x ind,box,lis]=make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps);


if dimx==1
    scalx=scal(1);
    [boxx1,lisx1,mxi]=make_box1(x(1,:),N,scalx,BOX1);
else
    [xx indx,boxx,lisx]=make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps);
end

if dimy==1
    scaly=scal(dimx);
    [boxy1,lisy1,myi]=make_box1(x(dimx+1,:),N,scaly,BOX1);
else
    [yy indy,boxy,lisy]=make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps);
end

dxy2=0;

for i=1:N
    xc=x(1:(dimx+dimy),ind(i));
    nn=neiK(x,dimx+dimy,0,dimx,ind(i),BOX,epsilon,K,box,lis,nn);
    epsx=0;
    for d=1:dimx
        for k=2:K+1
            dx=abs(xc(d)-x(d,nn(k)));
            if dx>epsx
                epsx=dx;
            end
        end
    end
    epsy=0;
    for d=dimx+1:dimx+dimy
        for k=2:K+1
            dy=abs(xc(d)-x(d,nn(k)));
            if dy>epsy
                epsy=dy;
            end
        end
    end
    
    if dimx>1
        nx2=neiE(xx,indx(i),0,dimx-1,dimx,BOX,epsilon,epsx,boxx,lisx);
    else
        nx2=neiE1(x(1,:),ind(i),scalx,BOX1,epsx,boxx1,lisx1,mxi);
    end
    if dimy>1
        ny2=neiE(yy,indy(i),0,dimy-1,dimy,BOX,epsilon,epsy,boxy,lisy);
    else
        ny2=neiE1(x(dimx+1,:),ind(i),scaly,BOX1,epsy,boxy1,lisy1,myi);
    end
    
    dxy2=dxy2+psi(nx2+1)+psi(ny2+1);
end
dxy2=dxy2/N;
mir=psi(N+1)+phi(K+1)-dxy2;

function ret=neiE(x, i, comp1, comp2, dim, bs, epsgrid, eps, box, lis)
comp1=comp1+1;
comp2=comp2+1;
ib=bs-1;
xx=x(1:dim,i);
ix=my_bitand(floor(xx(comp1)/epsgrid),ib);
iy=my_bitand(floor(xx(comp2)/epsgrid),ib);
jj=0; nx=0;
while eps>epsgrid*(jj-1)
    if jj==0
        step =1;
    else
        step = 2*jj;
    end
    for ix1=ix-jj:ix+jj
        ix2=my_bitand(ix1,ib);
        for iy1=iy-jj:step:iy+jj
            el=box(ix2+1,my_bitand(iy1,ib)+1);
            while el ~= -1
                dd=abs(xx(1)-x(1,el));
                for d=2:dim
                    dy=abs(xx(d)-x(d,el));
                    if dy>dd
                        dd=dy;
                        if dd>eps
                            break;
                        end
                    end
                end
                if dd<=eps
                    nx=nx+1;
                end
                el=lis(el);
            end
        end
    end
    for ix1=ix-jj:step:ix+jj
        ix2=my_bitand(ix1,ib);
        for iy1=iy-jj+1:iy+jj-1
            el=box(ix2+1,my_bitand(iy1,ib)+1);
            while el ~= -1
                dd=abs(xx(1)-x(1,el));
                for d=2:dim
                    dy=abs(xx(d)-x(d,el));
                    if dy>dd
                        dd=dy;
                        if dd>eps
                            break;
                        end
                    end
                end
                if dd<=eps
                    nx=nx+1;
                end
                el=lis(el);
            end
        end
    end
    jj=jj+1;
    if jj==(bs/2)
        break;
    end
end
if jj==(bs/2)
    for ix1=ix-jj:ix+jj
        ix2=my_bitand(ix1,ib);
        iy1=iy-jj;
        el=box(ix2+1,my_bitand(iy1,ib)+1);
        while el ~= -1
            dd=abs(xx(1)-x(1,el));
            for d=2:dim
                dy=abs(xx(d)-x(d,el));
                if dy>dd
                    dd=dy;
                    if dd>eps
                        break;
                    end
                end
            end
            if dd<=eps
                nx=nx+1;
            end
            el=lis(el);
        end
    end
    ix1=ix-jj;
    ix2=my_bitand(ix1,ib);
    for iy1=iy-jj+1:iy+jj-1
        el=box(ix2+1,my_bitand(iy1,ib)+1);
        while el ~= -1
            dd=abs(xx(1)-x(1,el));
            for d=2:dim
                dy=abs(xx(d)-x(d,el));
                if dy>dd
                    dd=dy;
                    if dd>eps
                        break;
                    end
                end
            end
            if dd<=eps
                nx=nx+1;
            end
            el=lis(el);
        end
    end
end
ret = nx-1;

function ret= neiE1(x, i, scal, bs, eps,  box, lis, mxi)
nx=0;
xc=x(i);
mp =floor((xc+eps)*scal);
if mp>bs
    mp=bs;
end
mm=floor((xc-eps)*scal);
if mm<0
    mm=0;
end
mi=box(mp+1);%+1
while mi>=0
    dd=x(mi)-xc;
    if abs(dd)<=eps
        nx=nx+1;
    end
    mi=lis(mi);%+1
end
if mm>=mp
    ret = nx-1;
    return;
end
mi=box(mm+1);%+1
while mi>=0
    dd=xc-x(mi);%+1
    if abs(dd)<=eps
        nx=nx+1;
    end
    mi=lis(mi);%+1
end

nx=nx+mxi(mp-1+1)-mxi(mm+1);%+1+1
ret =nx-1;

function nn=neiK(x, dim, comp1, comp2, i,bs, epsgrid, K, box, lis,nn)
comp1=comp1+1;
comp2=comp2+1;

ib=bs-1;

dn=zeros([1 K+1]);
xx=x(1:dim,i);

dn(1)=0;
dn(2:K+1)=1e30;

ix=my_bitand(floor(xx(comp1)/epsgrid),ib);
iy=my_bitand(floor(xx(comp2)/epsgrid),ib);

jj=0;
while dn(K+1)>epsgrid*(jj-1)
    if jj==0
        step =1;
    else
        step = 2*jj;
    end
    for ix1=ix-jj:ix+jj
        ix2=my_bitand(ix1,ib);
        for iy1=iy-jj:step:iy+jj
            el=box(ix2+1,my_bitand(iy1,ib)+1);
            while el ~= -1
                if el ~=i
                    dd=abs(xx(1)-x(1,el));
                    for d=2:dim
                        dy=abs(xx(d)-x(d,el));
                        if dy>dd
                            dd=dy;
                        end
                    end
                    if dd<dn(K+1)
                        k=K;
                        while dd<dn(k+1)
                            if k<K
                                dn(k+2)=dn(k+1);
                                nn(k+2)=nn(k+1);
                            end
                            k=k-1;
                        end
                        dn(k+2)=dd;
                        nn(k+2)=el;
                    end
                end
                el=lis(el);
            end
        end
    end
    for ix1=ix-jj:step:ix+jj
        ix2=my_bitand(ix1,ib);
        for iy1=iy-jj+1:iy+jj-1
            el=box(ix2+1,my_bitand(iy1,ib)+1);
            while el ~= -1
                if el ~=i
                    dd=abs(xx(1)-x(1,el));
                    for d=2:dim
                        dy=abs(xx(d)-x(d,el));
                        if dy>dd
                            dd=dy;
                        end
                    end
                    dy=abs(xx(2)-x(2,el));
                    if dy>dd
                        dd=dy;
                    end
                    if dd<dn(K+1)
                        k=K;
                        while dd<dn(k+1)
                            if k<K
                                dn(k+2)=dn(k+1);
                                nn(k+2)=nn(k+1);
                            end
                            k=k-1;
                        end
                        dn(k+2)=dd;
                        nn(k+2)=el;
                    end
                end
                el=lis(el);
            end
        end
    end
    jj=jj+1;
    if jj==bs/2
        break;
    end
end
if jj==(bs/2)
    for ix1=ix-jj:ix+jj
        ix2=my_bitand(ix1,ib);
        iy1=iy-jj;
        el=box(ix2+1,my_bitand(iy1,ib)+1);
        while el ~= -1
            if el~=i
                dd=abs(xx(1)-x(1,el));
                for d=2:dim
                    dy=abs(xx(d)-x(d,el));
                    if dy>dd
                        dd=dy;
                    end
                end
                if dd<dn(K+1)
                    k=K;
                    while dd<dn(k+1)
                        if k<K
                            dn(k+2)=dn(k+1);
                            nn(k+2)=nn(k+1);
                        end
                        k=k-1;
                    end
                    dn(k+2)=dd;
                    nn(k+2)=el;
                end
            end
            el=lis(el);
        end
    end
    ix1=ix-jj;
    ix2=my_bitand(ix1,ib);
    for iy1=iy-jj+1:iy+jj-1
        el=box(ix2+1,my_bitand(iy1,ib)+1);
        while el ~= -1
            if el~=i
                dd=abs(xx(1)-x(1,el));
                for d=2:dim
                    dy=abs(xx(d)-x(d,el));
                    if dy>dd
                        dd=dy;
                    end
                end
                if dd<dn(K+1)
                    k=K;
                    while dd<dn(k+1)
                        if k<K
                            dn(k+2)=dn(k+1);
                            nn(k+2)=nn(k+1);
                        end
                        k=k-1;
                    end
                    dn(k+2)=dd;
                    nn(k+2)=el;
                end
            end
            el=lis(el);
        end
    end
end

function [box, lis, mxi]=make_box1(x, N, scal, bs)
box(1:bs+1)=-1;
mxi(1:bs+1)=0;
for i=1:N 
    ix=floor(x(i)*scal);
    lis(i)=box(ix+1); 
    box(ix+1)=i; 
    mxi(ix+1)=mxi(ix+1)+1;
end
for i=2:bs+1
    mxi(i)=mxi(i)+mxi(i-1);
end

function [x ind,box,lis]=make_box2ind(x, dim, N, comp1, comp2, bs, inveps)
comp1=comp1+1;
comp2=comp2+1;
ib=bs-1;
xx=x;
box(1:bs,1:bs)=-1;

for i=1:N
    ix=my_bitand(floor(x(comp1,i)*inveps),ib)+1;
    iy=my_bitand(floor(x(comp2,i)*inveps),ib)+1;
    lis(i)=box(ix,iy);
    box(ix,iy)=i;
end

% ix=my_bitand(floor(x(comp1,1:N)*inveps),ib)+1;
% iy=my_bitand(floor(x(comp2,1:N)*inveps),ib)+1;
% lis=box(sub2ind(size(box),ix,iy));
% box(sub2ind(size(box),ix,iy))=1:N;
i=0;
for ix=1:bs
    for iy=1:bs
        ixy=box(ix,iy);
        while ixy >=0
            i=i+1;
            for d=1:dim
                x(d,i)=xx(d,ixy);
            end
            ind(ixy)=i;
            ixy=lis(ixy);
        end
        box(ix,iy)=-1;
    end
end

for i=1:N
    ix=my_bitand(floor(x(comp1,i)*inveps),ib)+1;
    iy=my_bitand(floor(x(comp2,i)*inveps),ib)+1;
    lis(i)=box(ix,iy);
    box(ix,iy)=i;
end

% ix=my_bitand(floor(x(comp1,1:N)*inveps),ib);
% iy=my_bitand(floor(x(comp2,1:N)*inveps),ib);
% lis(1:N)=box(sub2ind(size(box),ix,iy));
% box(sub2ind(size(box),ix,iy))=1:N;


function ret = my_bitand(a,b)
ret=typecast(bitand(typecast(int32(a),'uint32'),b),'int32');

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