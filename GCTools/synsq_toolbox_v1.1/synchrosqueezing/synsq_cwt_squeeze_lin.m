%
% function [Tx,fs] = synsq_cwt_squeeze(Wx, w, t, nv, opt)
%
% Calculates synchrosqueezed transform of f on a logarithmic
% scale.  Used internally by synsq_cwt_fw.
%
% Input:
%   Wx: wavelet transform of x
%   w: estimate of frequency at locations in Wx (see synsq_cwt_fw code)
%   t: time vector
%   nv: number of voices
%   opt: options struct, not currently used
%
% Output:
%   Tx: synchrosqueezed output
%   fs: associated frequencies
%
% Note the multiplicative correction term f in synsq_cwt_squeeze_mex (and in
% the matlab equivalent code commented out), required due to the fact that
% the squeezing integral of Eq. (2.7), in, [1], is taken w.r.t. dlog(a).
% This correction term needs to be included as a factor of Eq. (2.3), which
% we implement here.
%
% A more detailed explanation is available in Sec. III of [2].
% Specifically, this is an implementation of Sec. IIIC, Alg. 1.
% Note the constant multiplier log(2)/nv has been moved to the
% inverse of the normalization constant, as calculated in synsq_adm.m
%
% 1. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: a
% tool for empirical mode decomposition", 2010.
%
% 2. E. Brevdo, N.S. FuÄ?kar, G. Thakur, and H-T. Wu, "The
% Synchrosqueezing algorithm: a robust analysis tool for signals
% with time-varying spectrum," 2011.
%  
% 3. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed wavelet transforms: An empirical mode
% decomposition-like tool", 2011
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
%---------------------------------------------------------------------------------
function [Tx,fs,AvgM] = synsq_cwt_squeeze_lin(Wx, w, t, nv, opt)
    dt = t(2)-t(1);
    dT = t(end)-t(1);
    
    % Maximum measurable frequency of data
    %fM = 1/(4*dt); % wavelet limit - tested
    fM = 1/(2*dt); % standard
    % Minimum measurable frequency, due to wavelet
    fm = 1/dT; % really
    %fm = 1/(2*dT); % standard

    [na, N] = size(Wx);
    as = 2^(1/nv) .^ [1:1:na]'*2; %remove *2 to see beyond Fs/2
    
    %*** uniform frequencies ***%
    f1=1/as(1);
    flast=1/as(end);
    f=linspace(f1,flast,length(as));
    as=1./f';
    %*** uniform frequencies ***%
    
    das = [1; diff(as)];
    lfm = log2(fm); lfM = log2(fM);
    %fs = 2.^linspace(lfm, lfM, na);
    fs = linspace(fm, fM, na);
    %dfs = [fs(1) diff(fs)];

    % Harmonics of diff. frequencies but same
    % magnitude have same |Tx|
    dfs = diff(fs);
    dfs=[dfs(1) dfs];
    Dfs=repmat(dfs',[1 N]);
    

    if norm(Wx,'fro') < eps,
        Tx = zeros(size(Wx));
    else
        Tx = zeros(size(Wx));
        AvgM = zeros(size(Wx));
        for b=1:N
            for ai=1:length(as)
                if ai==1
                    dak = abs(as(ai) - as(ai+1));
                else
                    dak = abs(as(ai) - as(ai-1));
                end
                if (isfinite(w(ai, b)) && (w(ai,b)>0))
                    % Find w_l nearest to w(a_i,b)
                    %  2.^(lfm + k*(lfM-lfm)/na) ~= w(a_i,b)
                    %k = 1 + floor(na/(lfM-lfm)*(log2(w(ai,b))-lfm));
                    k = 1 + floor(na/(fM-fm)*(w(ai,b)-fm));
                    if isfinite(k) && k>0 && k<=na
                        % Tx(k,b) = Tx(k, b) + fs(k)/dfs(k) * Wx(ai, b) * as(ai)^(-1/2);
                        %Tx(k,b) = Tx(k, b) + Wx(ai, b) * as(ai)^(-1/2);
                        Tx(k,b) = Tx(k, b) + Wx(ai, b) * (as(ai)^(-3/2))*dak; %linear! [3] eq 2.4
                        AvgM(k,b) = AvgM(k, b) + 1;
                    end
                end
            end % for ai...
        end % for b
        %Tx = 1/nv * Tx;
        Tx=Tx./Dfs;
        AvgM(AvgM==0)=1;
    end
end

