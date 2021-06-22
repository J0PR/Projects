%[mres,freqs,p,sig_thresh]=causalMatrix(y,optStruct);
%[mres,freqs,p]=causalMatrixPartiallyConditioned(y(:,end-200:end),optStruct);
mres2=mres.*~eye(20);
mres2(mres2<0.26)=0;
figure,h=sanePColor(mres2);
colormap hot
shading faceted
set(gca,'YDir','reverse');
%%% thin colorbar..
c=colorbar;
x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
set(gca,'position',x1)
%%% thin colorbar..
set(gca,'YTick',1:20);
set(gca,'XTick',1:20);
set(h,'EdgeColor','w');
title('PDC')
xlabel('From')
ylabel('To')
set(findall(gcf,'type','text'),'fontSize',7)
set(gca,'FontSize',7)
xSize = 7; ySize = 6*xSize /8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
circularCausal(mres2, {'1' '2' '3' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '19' '20'}, [],  {sum(mres2,1),sum(mres2,2)',sum(mres2>0,1)+sum(mres2>0,2)'},  {'copper','autumn','summer'}, {'outflow','inflow','degree'}, [5,5,5,5], {'AR1','AR2','AR3','AR4'})
set(findall(gcf,'type','text'),'fontSize',7)
set(gca,'FontSize',7)
xSize = 7; ySize = 6*xSize /8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
set(gcf, 'Color', 'w');