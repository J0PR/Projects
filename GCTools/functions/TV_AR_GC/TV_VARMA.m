function [phi, theta, C]=TV_VARMA(IMFs_Matrix,p,q)
%Self_predict: IMFs from the same variable don't predict each others. Only IMFs from the
%same scale. Variables must have the same number of IMFs then!
%y : nPoints*nVars
%% INIT
if nargin<2
    p=2;
    q=1;
end
[r, c]=size(IMFs_Matrix);
if r<c
    IMFs_Matrix=IMFs_Matrix';
end
[nPoints,nIMFs]=size(IMFs_Matrix);

initial_length=0.1*nPoints; % % of data to include in initial data.
if initial_length<20
    initial_length=20;
end
initial_data=IMFs_Matrix(1:initial_length,:);
%% VARMA estimation
%strv=VARMA4IMF_est(initial_data,p,q,nIMFs);
strv=VARMA_est(initial_data,p,q);

%% TV-VARMA estimation
%ret=tvKalman_4IMF(strv.lastPhi(:,:,2:end),strv.lastTheta(:,:,2:end),IMFs_Matrix');
res=tvKalman(strv.lastPhi(:,:,2:end),strv.lastTheta(:,:,2:end),IMFs_Matrix');

phi=res.phi;
theta=res.theta;
C=res.ECov;
