% This function plots a 3D line (x,y,z) encoded with scalar color data (c)
% using the specified colormap (default=jet);
%
% SYNTAX: h=cline(x,y,z,c,colormap);
%
% DBE 09/03/02

function h=my_cline(x,y,z,scale,cmap,opt);

% if nargin==0  % Generate sample data...
%   x=linspace(-10,10,101);
%   y=2*x.^2+3;
%   z=sin(0.1*pi*x);
%   z=zeros(size(x));
%   c=exp(z);
%   c=1:length(z);
%   w=z-min(z)+1;
%   cmap='jet';
% elseif nargin<4
%   fprintf('Insufficient input arguments\n');
%   return;
% elseif nargin==4
%   cmap='jet';
% end

if isnan(sum(x))
    c=zeros(size(x));
    ramp=1;
    for i=1:length(c)
        if ~isnan(x(i))
            c(i)=ramp;
            ramp=ramp+1;
        else
            c(i)=1;
            ramp=1;
        end
    end
else
    c=1:length(x);
end
cmap=colormap(cmap);                      % Set colormap
yy=linspace(min(c),max(c),size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',c);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
% Lot line segment with appropriate color for each data pair...
cmScale=rgb2hsv(cm');
if strcmp(opt,'light')
    cmScale=hsv2rgb([cmScale(:,1) repmat(scale,[size(cmScale,1) 1]) ones(size(cmScale,1),1)]);
    for i=1:length(z)-1
        %         cmScale=rgb2hsv(cm(:,i)');
        %         cmScale=hsv2rgb([cmScale(1) scale 1]);
        %         h(i)=line([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'color',cmScale,'LineWidth',1);
        h(i)=line([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'color',cmScale(i,:),'LineWidth',1);
    end
elseif strcmp(opt,'dark')
    cmScale=hsv2rgb([cmScale(:,1) ones(size(cmScale,1),1) repmat(scale,[size(cmScale,1) 1])]);
    for i=1:length(z)-1
%         cmScale=rgb2hsv(cm(:,i)');
%         cmScale=hsv2rgb([cmScale(1) 1 scale]);
        h(i)=line([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'color',cmScale(i,:),'LineWidth',1);
    end
else
    cmScale=hsv2rgb([cmScale(:,1) repmat(scale,[size(cmScale,1) 1]) ones(size(cmScale,1),1)]);
    for i=1:length(z)-1
%         cmScale=rgb2hsv(cm(:,i)');
%         cmScale=hsv2rgb([cmScale(1) scale 1]);
        h(i)=line([x(i) x(i+1)],[y(i) y(i+1)],[z(i) z(i+1)],'color',cmScale(i,:),'LineWidth',1);
    end
end

return
