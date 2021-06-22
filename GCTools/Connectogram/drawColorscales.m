function hh=drawColorscales(hh,pos,m,M)
%pos = 'left' 'right' 'top' 'bot'
r=0.2;
N=length(hh.Rin);
R=hh.Rout(end);

%see if r needs to be adjusted
if N*2*r+N*2*m>2*M+2*R
    r=(2*M+2*R-N*2*m)/(2*N);
end
xx=ones(1,N)*(R+M+m+r);
yy=ones(1,N)*(R+M/2)-([0:(N-1)]*(2*r+2*m));
offset=0;
switch pos
    case 'right'
        offset=0;
    case 'top'
        offset=pi/2;
    case 'left'
        offset=pi;
    case 'bot'
        offset=3*pi/2;
end
[Ts,Rs] = cart2pol(xx,yy);
Ts=Ts+offset;

for i=1:N
    drawColorCircle(Ts(i),Rs(i),r,hh.colormap(:,:,i),hh.yy(i,:),hh.title{i});
end

function drawColorCircle(Ts,Rs,R,colorm,vals,title)
% colorm=[colorm(round(1:size(colorm,1)/20:size(colorm,1)),:);colorm(end,:)];
% vals=[vals(round(1:size(vals,1)/20:size(vals,1))) vals(end,:)];
%colorm=[linspace(colorm(1,1),colorm(end,1),20)',linspace(colorm(1,2),colorm(end,2),20)',linspace(colorm(1,3),colorm(end,3),20)'];
%vals=linspace(min(vals),max(vals),20);
x=ones(1,size(colorm,1));
x = x/sum(x);
theta0 = pi/2;
maxpts = 5;

[xtitle,ytitle] = pol2cart(pi/2,2*R);
[offx,offy] = pol2cart(Ts,Rs);

inside=0;
insertVals=[round(1:length(x)/4:length(x)) length(x)];
for i=1:length(x)
  n = max(1,ceil(maxpts*x(i)));
  r = [0;ones(n+1,1);0]*R;
  theta = theta0 + [0;x(i)*(0:n)'/n;0]*2*pi;
  if inside,
    [xtext,ytext] = pol2cart(theta0 + x(i)*pi,.5*R);
  else
    [xtext,ytext] = pol2cart(theta0 + x(i)*pi,1.2*R);
  end
  [xx,yy] = pol2cart(theta,r);
  
  theta0 = max(theta);
 
  if xtext < 0
      horizAlign = 'right';
  else
      horizAlign = 'left';
  end
  
  
  %sum the offsets
  xx=xx+offx;
  yy=yy+offy;

  
  xtext=xtext+offx;
  ytext=ytext+offy;
  
  if sum(insertVals==i)~=0%imprime label também
      patch(xx',yy',colorm(i,:),'linestyle','non')
      text(xtext,ytext,num2str(vals(i),'%0.1f'),'HorizontalAlignment',horizAlign,'FontSize',7);
  else
      patch(xx',yy',colorm(i,:),'linestyle','non')
  end
end
xtitle=xtitle+offx;
ytitle=ytitle+offy;
text(xtitle,ytitle,title,'FontSize',7,'HorizontalAlignment','center');