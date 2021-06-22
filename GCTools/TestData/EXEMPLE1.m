%load optStruct
%change its fields
load('TestData\AR_data_Biosig\optStruct.mat')

%load signal into variable y
load('TestData\AR_data_Biosig\AR_signal.mat');
%complete analysis: conditional or pairwise
[mres,freqs,p,sig_thresh]=causalMatrix(y,optStruct);
%partially conditioned analysis (for small number observations and large
%number of vars
%[mres,freqs,p]=causalMatrixPartiallyConditioned(y(:,end-200:end),optStruct);

%plotting
plotAdjM(mres,'test',{'1' '2' '3' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '19' '20'})

%conectogram 1
circularCausal(mres, {'1' '2' '3' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '19' '20'}, [],  {sum(mres,1),sum(mres,2)',sum(mres>0.2,1)+sum(mres>0.2,2)'},  {'copper','autumn','summer'}, {'outflow','inflow','degree'}, [5,5,5,5], {'AR1','AR2','AR3','AR4'})
%resize..
xSize = 8; ySize = 6;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])

%conectogram 2
circularCausal(mres, {'1' '2' '3' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '19' '20'}, [],  {sum(mres,1),sum(mres,2)',sum(mres>0.2,1)+sum(mres>0.2,2)'},  {'copper','autumn','summer'}, {'outflow','inflow','degree'}, [], [])
%resize..
xSize = 8; ySize = 6;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])