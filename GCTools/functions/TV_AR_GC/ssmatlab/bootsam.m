
function [Ystar, t1] = bootsam(y,X,Z,G,W,T,H,ins,i,beta,nboot)
%**************************************************************************
%**************************************************************************
%        This function applies the bootstrap algorithm by Stoffer and
%        Wall(1991) described below corresponding to the model:
%
%        y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%        alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%        collapse is applied to get rid of the diffuse component.
%
%        INPUTS:
%        y:     an (n x p) matrix of observations;
%        X    : an (n*p x nbeta) matrix containing the X_t matrices; 
%               a  (p x nbeta) if it is time invariant;
%               it can be []
%        Z    : an (n*p x nalpha) matrix containing the Z_t matrices;
%               a  (p x nalpha) matrix if it is time invariant
%        G    : an (n*p x nepsilon) matrix containing the G_t matrices;
%               a  (p x nepsilon) matrix if it is time invariant
%        W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%               an (nalpha x nbeta) matrix if it is time invariant;
%               it can be []
%        T    : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%               an (nalpha x nalpha) matrix if it time invariant
%        H    : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%               an (nalpha x nepsilon) if it is time invariant
%        ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%             state information, according to array i below
%        i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%               where
%               cc   = nalpha if c is not missing (0 if c missing)
%               cw0  = number of columns in W_0 (0 if W_0 missing)
%               ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%               cca1 = number of columns in A_1 (0 if A_1 missing)
%        beta : estimated regression vector
%        nboot: number of bootstraps
%
%        OUTPUT:
%        Ystar: an (n-t1+1) x p x nboot matrix of bootstrapped samples
%        t1   : time point of the collapse of the diffuse Kalman
%               filter

%*************************************************************************
%                   Bootstrap sample 
%
% Stoffer and Wall (1991): "Bootstrapping State-Space Models: Gaussian
% Maximum Likelihood Estimation and the Kalman Filter"
%
% Innovations form of the state space model written in the super matrix 
% form:
% Xi_t = A_t * Xi_{t-1} + B_t * beta + C_t * e_t,
% where:
% Xi_t = [ alpha_{t+1|t}; Y_t ], (nalpha + p) x 1  matrix
% e_t  = Sigma_t1(-1/2) * eps : standardized innovations, p x 1 matrix
% A_t = [ T_t  0 ; Z_t  0 ]
% B_t = [ W_t ; X_t  ]        
% C_t = [ K_t * Sigma_t^(1/2} 0 ;  Sigma_t^(1/2}  0 ]
% 
%
% Algorithm for obtaining bootstrap sample Ystar:
%
% Step 1: construct standardized innovations: e
% Step 2: bootstrap sample from the standardized innovations: estar
% Step 3: compute Xistar_t according to the model:
%         Xistar_t = A_t * Xistar_{t-1} + B_t * beta + C_t * estar_t
% Step 4: Ystar_t =  Xistar_t(nalpha+1:nalpha+p)
% 
%**************************************************************************
% Written by Martyna Marczak, 23.01.2012
% Department of Economics (520G) 
% University of Hohenheim
% Schloss, Museumsfluegel
% 70593 Stuttgart, Germany
% Phone: + 49 711 459 23823
% E-mail: marczak@uni-hohenheim.de
%**************************************************************************
%**************************************************************************


ii=i;


% Step 1: Obtain standardized innovations

% [~,~,recrs,~,srecr,t1,A1,~,KG]=scakfff(y,X,Z,G,W,T,H,ins,ii,beta);
[Xtf,Ptf,recrs,recr,srecr,t1,A1,P1,KG]=scakfff(y,X,Z,G,W,T,H,ins,ii,beta);  

% Obtain initial values for y_t in the Xi_t recursion
y1 = y(t1,:)'; 


% Step 2:

[le, p] = size(recrs);      
randidx = zeros(le,nboot);  % matrix of random indexes
Estar = zeros(le,p,nboot); % matrix of bootstrapped stand. innovations

%the following two lines are to generate always the same bootstrap samples
seed=25;
stream=RandStream('mt19937ar','Seed',seed);
     %Instead of the generator "shr3cong" (Shift-register generator
     %summed with linear congruential generator, of approximate period
     %2^64) used in MATLAB 7.1, MATLAB 7.7 now uses by default "mt19937ar" 
     %(Mersenne twister, of approximate period 2^19936 -1).
     %
     % The subsequent calls of rand should be of the type:
     % rand(stream,m,n). 
     % Other option:
     % seed=25;
     % stream=RandStream('mt19937ar','Seed',seed);
     % RandStream.setDefaultStream(stream);
     % rand(m,n);
     % Alternatively, use the command:
     % rng('default') to generate default random stream with seed = 0.
     
% Bootstrap from the standardized innovations
for j = 1:nboot
    randidxi = ceil(rand(stream,1,le)*le); 
    randidx(:,j) = randidxi;
    Estar(:,:,j) = recrs(randidxi,:);
end

% Step 3 and 4:

Ystar = zeros(le,p,nboot);
[mz, nalpha] = size(Z);
[mt, junk] = size(T);

if (mz == p && mt == nalpha)
    A = [ T  zeros(nalpha,p); Z zeros(p, p) ]; % Time invariant super matrix A
end

if isempty(beta)
   beta = 0;
   nbeta = 1;
else nbeta = length(beta);
end

if isempty(W)
   W = zeros(nalpha,nbeta);
end

if isempty(X)
   X = zeros(p,nbeta);
end

[mw, junk] = size(W);
[mx, junk] = size(X);

if (mx == p && mw == nalpha)
    B = [ W; X ];  % Time invariant super matrix B
end

Xi0 = [ A1; y1];   % Initializing Xi vector

for j = 1:nboot
    ystar = zeros(le,p);
    Xi=Xi0;
    for i = 1:le
        
        % Row indexes for the matrices
        ialpha = 1+nalpha*(i-1):nalpha*i;
        ip = 1+p*(i-1):p*i;
        
        % Partly or fully time-varying super matrix A
        if (mz > p || mt > nalpha)
            if mz == p
                A = [ T(ialpha,:) zeros(nalpha,p); Z zeros(p, p) ];
            elseif mt == nalpha
                A = [ T  zeros(nalpha,p); Z(ip,:) zeros(p, p) ];
            else
                A = [ T(ialpha,:)  zeros(nalpha,p); Z(ip,:) zeros(p, p) ];
            end
        end
        
        % Partly or fully time-varying super matrix B
        if (mx > p || mw > nalpha)
            if mx == p
                B = [ W(ialpha,:); X ];
            elseif mw == nalpha
                B = [  W;  X(ip,:) ];
            else
                B = [ W(ialpha,:);  X(ip,:)];
            end
        end
        
        %Compute Sigma_t^(1/2) and KGS_t = K_t * Sigma_t^(1/2}
        SigmaHlf = (chol(srecr(ip,:)))';
        KGS = KG(ialpha,:)*SigmaHlf;
        
        %Super matrix C
        C = [ KGS;  SigmaHlf ];
        
        estar = Estar(i,:,j)'; 

        Xi =  A*Xi + B*beta + C*estar;  
        
        ystar(i,:) = Xi(nalpha+1:end)';
    end
    Ystar(:,:,j) = ystar;   
    % NOTE: The bootstrapped sample is shorter than the original sample by
    % the number of observation up to the collapse of the diffuse Kalman
    % filter. Amendment of initial values can substantially distort 
    % the results, hence it is better not to extend the bootstrapped
    % sample.
end
end