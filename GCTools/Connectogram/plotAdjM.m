function plotAdjM(mres,title_str,label_cell)

figure,h=sanePColor(mres);
colormap hot
shading faceted
set(gca,'YDir','reverse');

set(gca,'YTick',1:size(mres,1));
set(gca,'XTick',1:size(mres,1));
set(gca,'YTickLabel',label_cell);
set(gca,'XTickLabel',label_cell);
xticklabel_rotate
%%% thin colorbar..
c=colorbar;
x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
set(gca,'position',x1)
%%% thin colorbar..
set(h,'EdgeColor','w');
title(title_str)
xlabel('From')
ylabel('To')
set(findall(gcf,'type','text'),'fontSize',7)
set(gca,'FontSize',7)
xSize = 7; ySize = 6*xSize /8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
set(gcf, 'Color', 'w');