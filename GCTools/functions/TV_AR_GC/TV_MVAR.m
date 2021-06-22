function [phi, C]=TV_MVAR(data,p)
%Self_predict: IMFs from the same variable don't predict each others. Only IMFs from the
%same scale. Variables must have the same number of IMFs then!
%y : nPoints*nVars
%% INIT
if nargin<2
    p=2;
end
[r, c]=size(data);
if r<c
    data=data';
end
%% TV-MVAR estimation
inp_model.data = data;
inp_model.order = p;
[phi, C] = DEKF3(inp_model);
