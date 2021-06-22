function circularCausal(M, varsLabels, mainCmap, valsList, cmapList, valsLabels, superGroupsLen, superGroupLabel)
% Circular Connectogram for adjacency matrices
% M: Nvars*Nvars adjacency matrix
% varsLabels: cell of length Nvars with strings for the names of each variable
% mainCmap: colormap for indicating directionality in the connection
%   curves. For custom colormaps, type their name without file extention. To
%   use black arrows instead of colored lines, leave this option empty.
% valsList: cell with Nvals arrays of length Nvars for additional values for each
%   variable. Ex: network metrics like inflow or outflow density,
%   clustering coeff, etc..
% cmapList: cell with Nvals strings with the colormap for the scales of
%   each value present in valsList
% valsLabels: cell with Nvals strings with the names of each vale present
%   in valsList
% superGroupsLen: array with length for possible subgroups of variables.
%   CAUTION: sum(superGroupsLen) == Nvars.
% superGroupLabel: cell with labels for the previous subgroups
% Examples:
% using arrows...
% circularCausal(magic(5), {'a' 'b' 'c' 'd' 'e'}, [],  {[1:5],[11:15]},  {'summer','copper'}, {'metrics1','metrics2'}, [2,3], {'group1','group2'})
% using colored (with the 'cool' colormap) lines for direction
% circularCausal(magic(5), {'a' 'b' 'c' 'd' 'e'}, 'cool',  {[1:5],[11:15]},  {'summer','copper'}, {'metrics1','metrics2'}, [2,3], {'group1','group2'})
% 
% save figure using the export_fig package by: export_fig figure_name.pdf
% (or other extention..)
%
% Author: João Rodrigues (based in the Oleg Komarov's schemball function)


%text color
tcolor = 'k';
% Number of color shades/buckets (large N simply creates many perceptually indifferent color shades)
N      = 20;
% Points in [0, 1] for bezier curves: leave space at the extremes to detach a bit the nodes. 
% Smaller step will use more points to plot the curves.
t      = (0.025: 0.05 :1)';
if isempty(mainCmap)
    use_arrows=1;
else
    use_arrows=0;
    if ischar(mainCmap)
        if exist([mainCmap '.mat']) == 2 %custom coplormap saved in the root folder..
            tempStruct=load([mainCmap '.mat']);
            structNames=fieldnames(tempStruct);
            mainCmap=tempStruct.(structNames{1});
        end
    end
end
vals=cell2mat(valsList');%dim: numValues*NumVariables
if ~isempty(superGroupsLen)%separate subgroups with NaN cols and rows.
    if isempty(superGroupLabel)
        superGroupLabel=repmat({''},size(superGroupsLen));
    end
    startPos=1;
    j=1;
    for i=1:length(superGroupsLen)-1
        len=superGroupsLen(i);
        M=[M(:,1:(startPos+len-1)) NaN*ones(size(M,1),1) M(:,(startPos+len):end)];
        M=[M(1:(startPos+len-1),:);NaN*ones(1,size(M,2));M((startPos+len):end,:)];
        vals=[vals(:,1:(startPos+len-1)) NaN*ones(size(vals,1),1) vals(:,(startPos+len):end)];
        varsLabels={varsLabels{1:(startPos+len-1)} '' varsLabels{(startPos+len):end}};
        startPos=startPos+len+1;
        
        %%%%%%%%%insert spaces in supergroups%%%%%%%%%%
        superGroupLabel={superGroupLabel{1:j} '' superGroupLabel{(j+1):end}};
        j=j+2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
     M(end+1,:)=NaN;
     M(:,end+1)=NaN;
     vals(:,end+1)=NaN;
     varsLabels{end+1}='';
     superGroupLabel{end+1}='';
end

% Create figure
figure('renderer','zbuffer','visible','on')%visible off
axes('NextPlot','add')
set(gca,'Visible','off')

%% Index causal relations..
tf=M>0;
sz = size(M);

% Index correlations into bucketed colormap to determine plotting order (brightest to darkest)
[n, isrt] = histc(M(tf), linspace(min(M(tf)),max(M(tf)) + eps(100),N + 1));

% Retrieve pairings of nodes
[row,col] = find(tf);

tau   = 2*pi;

% Positions of nodes on the circle starting from (0,-1), for label orientation
step  = tau/sz(1);
theta = -.25*tau : step : .75*tau - step;

% Get cartesian x-y coordinates of the nodes
x     = cos(theta);
y     = sin(theta);

%% PLOT BEZIER CURVES 
% Calculate Bx and By positions of quadratic Bezier curves with P1 at (0,0)
% B(t) = (1-t)^2*P0 + t^2*P2 where t is a vector of points in [0, 1] and determines, i.e.
% how many points are used for each curve, and P0-P2 is the node pair with (x,y) coordinates.

t2  = [1-t, t].^2;
s.l = NaN(N,1);

%% LOOP per color bucket
for pos = 1:N
    %pos = plotorder(c);
    idx = isrt == pos;
    
    if nnz(idx)
        Bx     = [t2*[x(col(idx)); x(row(idx))]; NaN(1,n(pos))];
        By     = [t2*[y(col(idx)); y(row(idx))]; NaN(1,n(pos))];
        if use_arrows
            colorVal=repmat(1-pos/N,[1 3]);%1-pos/N para ir do mais branco para o mais escuro (em fundo branco). Caso o fundo seja escuro usar pos/N.
            plot(Bx(:),By(:),'Color',colorVal,'LineWidth',1);
            startArrows=[Bx(end-2,:)' By(end-2,:)'];
            endArrows=[Bx(end-1,:)' By(end-1,:)'];
            arrow(startArrows,endArrows,'EdgeColor',colorVal,'FaceColor',colorVal)
        else
            my_cline(Bx(:),By(:),zeros(size(By(:))),pos/N,mainCmap,'light');
        end
    end
end

%% draw rings
s.rings=[];
for i=1:size(vals,1)
    s.rings=drawRing(theta,vals(i,:),i,0,0.2,cmapList{i},valsLabels{i},s.rings);
end

%% draw colorCircles
drawColorscales(s.rings,'right',0.2,0.4);

%% PLACE TEXT LABELS to read 'left to right'

Rlabels=s.rings.Rout(end);
[xlabel,ylabel]=pol2cart(theta,Rlabels);
ipos       = x > 0;
s.t        = zeros(sz(1),1);
s.t( ipos) = text(xlabel( ipos)*1.01, ylabel( ipos)*1.01, varsLabels( ipos),'Color',tcolor);
set(s.t( ipos),{'Rotation'}, num2cell(theta(ipos)'/tau*360))
s.t(~ipos) = text(xlabel(~ipos)*1.01, ylabel(~ipos)*1.01, varsLabels(~ipos),'Color',tcolor);
set(s.t(~ipos),{'Rotation'}, num2cell(theta(~ipos)'/tau*360 - 180),'Horiz','right')

% maxL=0;
% for i=1:length(varsLabels)
%     tmpLbl=varsLabels{i};
%     tmpL=length(tmpLbl);
%     if tmpL>maxL
%         maxL=tmpL;
%     end
% end

for i=1:length(varsLabels)
    tmpLbl=varsLabels{i};
    tmpL(i)=length(tmpLbl);
end
maxL=median(tmpL);

s.Text.Rin=Rlabels;
s.Text.Rout=Rlabels+maxL*0.1;

if ~isempty(superGroupsLen) && ~isempty(superGroupLabel)
    %% subgroups of variables
    % define angles for variables that ar sets of these subvariables..
    j=1;
    aux=0;
    sumAux=0;
    for i=1:size(vals,2)
        if(isnan(vals(1,i)))
            if i~=1
                stheta(j)=sumAux/aux;
                aux=0;
                sumAux=0;
                j=j+1;
            end
            stheta(j)=theta(i);
            j=j+1;
        else
            if i==size(vals,2)
                stheta(j)=sumAux/aux;
                aux=0;
                sumAux=0;
                j=j+1;
            else
                sumAux=sumAux+theta(i);
                aux=aux+1;
            end
        end
    end
    sx     = cos(stheta);
    sy     = sin(stheta);
    
    
    % PLACE TEXT LABELS for subgroups
    Rlabels=s.Text.Rout;
    [xlabel,ylabel]=pol2cart(stheta,Rlabels);
    ipos       = sy > 0;
    s.t        = zeros(sz(1),1);
    s.t( ipos) = text(xlabel( ipos)*1.2, ylabel( ipos)*1.2, superGroupLabel( ipos),'Color',tcolor,'FontWeight','bold','HorizontalAlignment','center');
    set(s.t( ipos),{'Rotation'}, num2cell(stheta(ipos)'/tau*360 - 90))
    s.t(~ipos) = text(xlabel(~ipos)*1.2, ylabel(~ipos)*1.2, superGroupLabel(~ipos),'Color',tcolor,'FontWeight','bold','HorizontalAlignment','center');
    %set(s.t(~ipos),{'Rotation'}, num2cell(stheta(~ipos)'/tau*360 + 90),'Horiz','right')
    set(s.t(~ipos),{'Rotation'}, num2cell(stheta(~ipos)'/tau*360 + 90))
    
    maxL=0;
    for i=1:length(superGroupLabel)
        tmpLbl=superGroupLabel{i};
        tmpL=length(tmpLbl);
        if tmpL>maxL
            maxL=tmpL;
        end
    end
    s.sText.Rin=Rlabels;
    s.sText.Rout=Rlabels+maxL*0.1;
end
%% background color and size
set(gcf, 'Color', 'w');

%xSize = 8; ySize = 5.4;
%xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
%set(gcf,'Position',[0 0 xSize*50 ySize*50])


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