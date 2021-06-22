function plotAdjM_freqs(mres,f,title_str,label_cell)
% mres 4D adj matrix
% t: time in seconds
% f: frequency Hz
% label_cell cell with label strings
figure,
plot_idx=1;
nvars=size(mres,1);
mres=mres.*repmat(~eye(size(mres,1)),[1 1 size(mres,3)]);
max_val=max(mres(:));
%mres(mres<0.2)=0;
for i=1:nvars
    for j=1:nvars
    
            ax(plot_idx)=subplot(nvars,nvars,plot_idx);
            area(f,abs(squeeze(mres(i,j,:))))
            ylim([0 max_val]);
            set(gca,'FontSize',7)
            if j==1 % col 1
                ylabel({label_cell{i}});
            end
            if i==nvars %last line
                xlabel({'f (Hz)',label_cell{j}});
            end
      

        plot_idx=plot_idx+1;
    end
end

suplabel('From');
suplabel('To','y');

set(findall(gcf,'type','text'),'fontSize',7)
set(gca,'FontSize',7)

suptitle(title_str);


xSize = 13; ySize = 6*xSize /8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
set(gcf, 'Color', 'w');