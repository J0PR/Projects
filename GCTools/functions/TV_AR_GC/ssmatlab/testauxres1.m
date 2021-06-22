function testauxres1(fid,tname,DS,MDS,str,result,dist,cw,gflag,out)
%**********************************************************************
%**********************************************************************
%   This function plots the smoothed disturbances of a univariate
%   structural time series model, which are specified by the user
%
%   INPUTS:
%   fid  : identifier of the file for writing the test results
%   tname: name of the underlying time series
%   p    : number of variables in the (n x p) observation matrix y
%   DS   : an n x (p+nalpha) matrix containing the estimated 
%          [G_t;H_t]*eps_{t|n}
%   MDS  : an (n*(p+nalpha)) x (p+nalpha) matrix containing the Mse of DS
%   str  : structure containing model information
%   dist : logical 1 x 3 array index,
%          used to specify which smoothed should be plotted:
%          1st element: 1 - plot irregular
%                       0 - do not plot irregular;
%          2nd element: 1 - plot trend disturbance
%                       0 - do not plot trend disturbance;
%          3rd element: 1 - plot slope disturbance
%                       0 - do not plot slope disturbance
%          Alternatively, in the case of a univariate structural time
%          series model selected disturbances can specified by one of the 
%          strings given below, or a cell array,if more disturbances are
%          to be plotted
%          'irreg'- irregular
%          'level'- trend disturbance
%          'slope'- slope slope disturbance
%          Example:
%          dist = 'level'
%          dist = {'level','irreg'}
%  cw    : critical point of the t distrubution
%  gflag : 1, pause between displaying figures
%          0, no pause
%  out   : write values for auxiliary residuals larger than out
%*************************************************************************
%**************************************************************************
% Written by Martyna Marczak,  26.01.2012
% Department of Economics (520G) 
% University of Hohenheim
% Schloss, Museumsfluegel
% 70593 Stuttgart, Germany
% Phone: + 49 711 459 23823
% E-mail: marczak@uni-hohenheim.de
%**************************************************************************
%**************************************************************************


% Checks for inconsistencies

if ~islogical(dist) && ~ischar(dist)  && ~iscell(dist)
    error('dist must be either a logical vector or a string/cell array')
end

% Transform a string to a cell array
if ischar(dist)
   dist = {dist};
end

if iscell(dist)
    for i = 1:length(dist)
       if ~any(strcmpi(dist{i},{'irreg','level','slope'}))
           error([dist{i} ' is not an allowed string'])
       end
    end
end

% Check for the existence of a univariate structural time series model
% by checking the existence of trend only

if any(strcmpi('trend',fieldnames(str)))
   ustsm = 1;
else 
   ustsm = 0;
end

if ustsm == 0
    error('The model is not a univariate structural time series model. Use logical array to indicate disturbances')
elseif ustsm == 1
% Transform the cell array into an equivalent logical vector
        if iscell(dist)
            d = zeros(1,3);
            if any(strcmpi('irreg',dist))
                d(1) = 1;
            else d(1) = 0;
            end
            if any(strcmpi('level',dist))
                d(2) = 1;
            else d(2) = 0;
            end
            if any(strcmpi('slope',dist))
                d(3) = 1;
            else d(3) = 0;
            end
            dist = logical(d);
        end
    
% Check whether the disturbances exist in the model under consideration or
% if they have standard deviation of zero (estimated or fixed)

        stord = str.stord;
        conc = str.conc;
        xvf = result.xvf;
        xf = result.xf;
        pfix = str.pfix;
        pvar = str.pvar;
        
        % Create vector with all parameters after estimation
        lx = length(xvf) + length(xf);
        x = zeros(lx,1);
        x(pfix) = xf;
        x(pvar) = xvf;
        
        sd = zeros(1,3);
           
        for i = 1:3
            if any(stord == i) 
                j = (stord == i);
                if x(j) == 0    % the estimated or fixed variance is equal to zero
                   sd(i) = 0;
                else sd(i) = 1;
                end
            elseif (conc ~= i) && ~any(stord == i) % the corresponding component is constant
                sd(i) = 0;
            else sd(i) = 1;
            end 
        end
        
        if (dist(1) == 1 && sd(1) == 0)
            error('Smoothed irregular is equal to zero')   
        end
        if (dist(2) == 1 && sd(2) == 0)
            error('Smoothed level disturbance is equal to zero')    
        end
        if (dist(3) == 1 && sd(3) == 0)
           error('Smoothed slope disturbance is equal to zero')    
        end
end



datei = str.datei;
freq = datei.freq;
[junk,nalpha]=size(str.T);
[n ,eps] = size(DS);



p = 1;      % number of variables in the observation matrix
ne = p+nalpha;
ndist = sum(dist);
SE = zeros(n,ne);

% Compute standardized smoothed disturbances
for i = 1:n
    ine = (i-1)*ne+1:i*ne;
    
    DDS = DS(i,:);
    MDDS = MDS(ine,:);
    % Standardizing
    SE(i,:) = DDS./(sqrt(diag(MDDS)))';
end

% Extend the selection vector dist by appending a zero vector to obtain a
% selection vector of the same length as the number of elements of the
% vector [G;H]*epsilon, i.e. p+nalpha = 1+nalpha

 dd = false(1,p+nalpha-3);
 dist = [dist dd];
 SE = SE(:,dist);  
 
 
% Assign names to the selected disturbances ordered as follows:
% irregular, level, slope

 dname = {'Irregular','Level','Slope'};
 
 k = dist(1:3) == 1;
 dname = dname(k(:));
 

    
% Create vector of times

beg_yr = datei.beg_yr;
beg_per = datei.beg_per;
beg_t = (beg_yr+(beg_per/freq));
end_t = (beg_yr+(beg_per/freq)) + (n/freq) - (1/freq);
t = (beg_t:(1/freq):end_t)';

% Show fractions as Q1 to Q4 for quarterly data or M1 to M12 for monthly
% data

tf = cell(freq,1);

if freq ~= 1
   tf = cell(freq,1);
   for i = 1:freq
       if freq == 4
          tf{i} = ['.Q',num2str(i)];   
       elseif freq == 12
          tf{i} = ['.M',num2str(i)];
       end
   end
end

times = cell(n,1);

for i = 1:n
    if freq == 1
        times{i} = num2str(t(i));
    elseif freq == 4 || freq == 12  
          q = t(i)-floor(t(i));
          % Modified 12.10.2012
          if q == 0
              if freq == 4
                  per = 4;
              else
                  per = 12;
              end
          else
             per = round(q*freq);
          end
          if per == 0
             per = freq;
          end
          if per ~= freq
             times{i} = [num2str(floor(t(i))),tf{per}];
          else
             times{i} = [num2str(floor(t(i))-1),tf{per}];
          end
    end
end

%*************************************************************************
%                Plot auxiliary residuals

vd = zeros(1,ndist);
std = zeros(1,ndist);

 for i = 1:ndist
     
    vd(1,i) = var(SE(:,i),1);   
    std(1,i) = sqrt(vd(1,i));
    
    figure
    bd=ones(n,1)*(cw*std(1,i));                          % confidence bands
    plot(t,SE(:,i),t,zeros(n,1),t,bd,'r',t,-bd,'r')      % plot residuals
   
    title('Auxiliary resididuals')
    legend(dname{i})
    
    if gflag == 1
    disp('strike any key when ready')
    pause
    end
 end
 
 %*************************************************************************
 %                            Tests
 %
 % Harvey A.C. and Koopman S.J. (1992) "Diagnostic Checking of Unobserved-Components
 % Time Series Model". Journal of Business & Economic
 % Statistics, 10, pp.377-389
 %
 % Correction factor: eq.(3.2)
  lag = max(sqrt(n),20);    % number of autocorrelations, truncation limit
                            % proposed in the STAMP manual
  corrf2 = zeros(1,ndist);  % correction factor for alpha=2
  corrf3 = zeros(1,ndist);  % correction factor for alpha=3
  corrf4 = zeros(1,ndist);  % correction factor for alpha=4
  K = zeros(ndist,2);       % Excess kurtosis statistic and its p-value
  N = zeros(ndist,2);       % Bowman-Shenton normality statistic and its p-value
  
  
  for i = 1:ndist
      m2 = var(SE(:,i),1);                % sample variance
      m3 = moment(SE(:,i),3);             % sample 3th moment
      m4 = moment(SE(:,i),4);             % sample 4th moment
      [c0,cv,r] = autcov(SE(:,i),lag,1);  % autocorrelations
      corrf2(i) = sum(2*(r.^2));
      corrf3(i) = sum(2*(r.^3));
      corrf4(i) = sum(2*(r.^4));
      kk = ((m4/m2^2)-3)/(sqrt(24*corrf4(i)/n)); % excess kurtosis statistic 
      K(i,1) = kk;
      pv = 1-normcdf(kk,0,1);                 % p-value of kk
      K(i,2) = pv;
      % Bowman-Shenton normality test statistic
      nn = (n*(m3/(m2^(3/2)))^2)/(6*corrf3(i)) + ...
           (n*((m4/m2^2)-3)^2)/(24*corrf4(i));
      N(i,1) = nn;
      pv = 1-chi2cdf(nn,2);                 % p-value of nn
      N(i,2) = pv;
  end
 
  
 %*************************************************************************
 %                     Print test results
 
 
 fprintf(fid,'\nSeries: %s\n', tname);
 fprintf(fid,'\n==================================================\n');
 for i = 1:ndist
     fprintf(fid,'           Auxiliary residuals: %s\n', dname{i});
     fprintf(fid,'\nAutocorrelation correction factor\n');
     clear in
     z = [ corrf2(i); corrf2(i); corrf3(i) ];
     in.rnames = char('    ','kappa(2)','kappa(3)','kappa(4)');
     in.fmt = char('%12.4f');
     in.fid = fid;
     mprint(z,in);
     fprintf(fid,'\n');
     
     clear in
     z = [K(i,1)  K(i,2);N(i,1) N(i,2)];
     in.cnames = char('   ','    P-values');
     in.rnames = char('   ','K: excess kurtosis test', 'N: Bowman-Shenton normality test');
     in.fmt = char('%12.4f','%12.4f');
     in.fid = fid;
     mprint(z,in)
     fprintf(fid,'\n');
     
     if any(abs(SE(:,i)) > out)
         j = abs(SE(:,i)) > out; 
         se = SE(:,i);  
         ose = se(j);
         date = {times{j}}'; 
        fprintf(fid,'   Date  '); fprintf(fid,['    |Outlier value|> ', num2str(out)]);
       for l = 1:length(ose)
           fprintf(fid,'\n%s    %15.2f', date{l}, ose(l));
       end
       
     end
     
     fprintf(fid,'\n--------------------------------------------------\n');
 end
 
end


