function plotAdjM_alt(mres,title_str,label_cell)

figure,h=imagesc(mres);
%%% thin colorbar..
c=colorbar;
x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
set(gca,'position',x1)
%%% thin colorbar..
set(gca,'YTick',1:size(mres,1));
set(gca,'XTick',1:size(mres,1));
set(gca,'YTickLabel',label_cell);
set(gca,'XTickLabel',label_cell);
title(title_str)
xlabel('From')
ylabel('To')

