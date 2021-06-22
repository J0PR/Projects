function hh = drawRing(phi,x,ring_num,Roffset,thickness,cmap,title,hh)

R1=0.96;
Rin=R1;
defaultThickness=0.2;
if ~isempty(hh)&&ring_num~=1
    Rin=hh.Rout(ring_num-1);
elseif ring_num~=1
    Rin=R1+(ring_num-1)*defaultThickness;
end
Rin=Rin+Roffset;
Rout=Rin+thickness;
hh.Rin(ring_num)=Rin;
hh.Rout(ring_num)=Rout;

% xs=x-min(x);
% xs=xs/max(xs);

cmap=colormap(cmap);                      % Set colormap
yy=linspace(min(x),max(x),size(cmap,1));  % Generate range of color indices that map to cmap
xm = spline(yy,cmap',x);                  % Find interpolated colorvalue
xm(xm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
xm(xm<0)=0;

hh.yy(ring_num,:)=yy;

maxpts = 100;
phi=phi+2*pi;
half_way=abs((phi(2)-phi(1))/2);
%fronteiras de phi(i) sao de phi(i)-half_way a phi(i)+half_way
for i=1:length(phi)
    n=max(1,ceil(maxpts/length(phi)));
    R = [ones(1,n)]*Rout;
    r = [ones(1,n)]*Rin;
    theta=linspace(phi(i)-half_way, phi(i)+half_way,n);
    [xx,yy] = pol2cart(theta,r);
    [XX,YY] = pol2cart(theta,R);
    if ~sum(isnan(xm(:,i)'))
        P = patch([XX flipdim(xx,2)],[YY flipdim(yy,2)],xm(:,i)');
    end
end
%guardar os dados das escalas de cor para desenhar as colorcicles no fim
hh.colormap(:,:,ring_num)=cmap;
hh.title{ring_num}=title;

