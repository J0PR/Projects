
function Cc = dispcomp(KKP,str,comp,varargin)

%************************************************************************
%************************************************************************
% This function displays and optionally plots selected components 
% of a univariate unobserved components model corresponding either to 
% the structural time series (STS) approach or to the ARIMA model based
% (AMB) approach
%
%    INPUTS:
%    REQUIRED inputs:
%    KKP   : smoothed state vector
%    str   : structure of model based on STS or AMB approach
%    comp  : logical array index, used to specify which elements of the 
%            smoothed states are to be output
%            Alternatively, in the case of a univariate unobserved
%            components model selected component can specified by one of the 
%            strings given below, or a cell array, if more components
%            should be displayed.
%            The order of components in the cell array does not play any role.
%            Allowed strings for structural time series (STS) model:
%               'level'- trend component
%               'slope'- slope component
%               'seas' - seasonal component
%               'cycle'- cycle component
%               'ar'   - AR component
%            Examples:
%               comp = 'level'
%               comp = {'level','cycle'}
%            Allowed strings corresponding to the ARIMA model based (AMB)
%            approach:
%               'trendcycle'- trend-cycle component
%                    'trend'- trend component
%                    'cycle'- cycle component
%                    'seas' - seasonal component 
%                    'tran' - transitory component
%            Examples:
%               comp = 'trendcycle'
%               comp = {'trend','seas'}
%  
%     OPTIONAL inputs (if components are to be plotted):
%     datei  : calendar structure (required input for plots)
%       ind  : scalar or a array with scalars specifying which of 
%              the components selected by comp should be plotted;
%              if ind is not input to dispcomp, all components given by
%              comp are plotted
%     vnames : string or array with strings to label components to be plotted;
%              number of names must be equal to length of comp, if ind is
%              not input to dispcomp;
%              number of names must be equal to length of ind, if ind is
%              input to dispcomp;
%              if vnames is not input to dispcomp, strings given in comp a
%              are used to label components to be plotted
%
%     OUTPUT:
%       Cc   : (n x ncomp) matrix with components that should be displayed
%     
%       If more than three arguments are passed to dispcomp, the function 
%       additionally plots the components 
% 
%     Examples:
%         Cc = dispcomp(KKP,str,'trend','Trend component')
%         Cc = dispcomp(KKP,str,{'trend','cycle'},2)
%         CC = dispcomp(KKP,str,{'trend','cycle'},2,'Cyclical component')
%
%**************************************************************************
% Written by Martyna Marczak, 23.01.2012
% Department of Economics (520G) 
% University of Hohenheim
% Schloss, Museumsfluegel
% 70593 Stuttgart, Germany
% Phone: + 49 711 459 23823
% E-mail: marczak@uni-hohenheim.de
%*************************************************************************
%*************************************************************************


% Check for inconsistencies
if ~islogical(comp) && ~ischar(comp)  && ~iscell(comp)
    error('comp must be either a logical vector or a string/cell array')
end


% Transform a string to a cell array
if ischar(comp)
   comp = {comp};
end

% Check for the existence of a univariate unobserved components model
% by checking the existence of trend only

if iscell(comp)
    if (~any(strcmpi('trend',fieldnames(str)))) && (~any(strcmpi('compcd',fieldnames(str))))
       error('The model is not a univariate unobserved components model. Use logical array to indicate components')
    end
end

if isfield(str,'stra')
    freq = str.stra.freq;
else
    freq = str.freq;
end

Z = str.Z;

% Check the existence of components in the model under consideration

if iscell(comp)
  
   if isfield(str,'trend')
     
      % STS approach
      if (any(strcmpi('level',comp)) && (str.trend==0 || str.trend==-1))
          error('Level is not included in the state vector')    
      end

      if (any(strcmpi('slope',comp)) && (str.slope==0 || str.slope==-1))
         error('Slope is not included in the state vector')    
      end

      if (any(strcmpi('seas',comp)) && str.seas==0)
         error('Seasonal is not included in the state vector')   
      end

      if (any(strcmpi('cycle',comp)) && str.cycle==0)
         error('Cycle is not included in the state vector')    
      end

      if (any(strcmpi('ar',comp)) && str.arp==0)
         error('AR is not included in the state vector')   
      end
      
      if str.trend == 1         % stochastic level
         level = 1;
      elseif str.trend == 2     % Butterworth tangent
         nlev = 2; 
         level = ones(1,nlev);
      else level = [];          % no or constant level
      end

      if str.slope == 1         % stochastic slope
         slope = 2;
      else slope = [];          % no or constant slope
      end

      if str.seas ~= 0
         if str.seas == 4
            nseas = freq;
         else nseas = freq-1;   % Butterworth tangent seasonality
         end
         seas = repmat(3,1,nseas);
      else seas = [];           % no seasonals
      end

      if str.cycle ~= 0
         cycle = [4 4];
      else cycle = [];          % no cycle
      end


      if str.arp ~= 0
         nar = str.arp;
         ar = repmat(5,1,nar); 
      else ar = [];             % no AR component
      end

      compall = [level slope seas cycle ar];

   else
       
      % AMB approach
      
      if any(strcmpi('trendcycle',comp))
          if  ~isfield(str,'compcd') || isfield(str,'compf')
             error('Trend-cycle is not included in the state vector')  
          end
      end
      
      if any(strcmpi('trend',comp)) 
          if ~isfield(str,'compcd') && ~isfield(str,'compf')
              error('Trend is not included in the state vector') 
          end
      end
      
      if any(strcmpi('cycle',comp))
          if ~isfield(str,'compcd') && ~isfield(str,'compf')
              error('Cycle is not included in the state vector')    
          end
      end
      
      if any(strcmpi('seas',comp)) 
          if ~isfield(str,'compcd')
              error('Seasonal is not included in the state vector')   
          end
      end
      
      if any(strcmpi('tran',comp)) 
          if ~isfield(str,'compcd') || (isfield(str,'compcd') && isempty(str.compcd.rt))
             error('Transitory component is not included in the state vector')
          end
      end
      
      if  isfield(str,'compcd') && ~isfield(str,'compf')
          ptden = str.compcd.ptden;
          np = length(ptden);
          trendcycle = ones(1,np);
      else
          trendcycle = [];
      end
      
      
      if isfield(str,'compcd')
          stden = str.compcd.stden; 
          ns = length(stden);
          seas = ones(1,ns)*4;
      else
          seas = []; ns=0;
      end
      
      
      if isfield(str,'compcd') && ~isempty(str.compcd.rt)
          rt = str.compcd.rt;
          nr = length(rt);
          tran = ones(1,nr)*5;
      else
          tran = []; nr=0;
      end
      
       
      if isfield(str,'compcd') && isfield(str,'compf')
          ptden = str.compcd.ptden;
          den = str.compf.den;
          nsp = length(den) + length(ptden)-1;
          trend = ones(1,nsp)*2;
          nc = length(Z) - (nsp+ns+nr);
          cycle = ones(1,nc)*3;
      else
          trend = []; cycle = [];
      end
      
     
      compall = [trendcycle trend cycle seas tran];
      
   end

   lcomp = length(comp);
   c = zeros(lcomp,length(compall));
 
   for i = 1:lcomp
       d = eval(comp{i});
       k = find(d(1) == compall);
       l = Z(k);
       c(i,k:k+length(d)-1) = l;
   end
   c = logical(c);
   
end

[lkkp,junk] = size(KKP);

Cc = zeros(lkkp,lcomp);

for i = 1:lcomp
    ci = c(i,:);
    com = sum(KKP(:,ci),2);
    Cc(:,i) = com;
end



% Check length of the variable input list
nargs = length(varargin);

% If there are more than three inputs to dispcomp, plot selected components

if nargs > 0
    
   datei = varargin{1};   % calendar structure as required input
   if ~isstruct(datei)
      error('calendar structure is required as input')
   end
   
   if nargs == 1
      for i = 1:lcomp
          figure 
          tsplot(Cc(:,i),datei,comp(i))
          pause
      end
      
   elseif nargs > 1
       
          if isscalar(varargin{2}) || ischar(varargin{2})
              varargin{2} = {varargin{2}};
          end
          
          if ~iscellstr(varargin{2})
              
              indc = varargin{2};
              for i=1:length(indc)
                  if ~isscalar(indc{i})
                      error('array must contain either only strings or only scalars')
                  end
              end
              
              indc = cell2mat(indc);
              if any(indc > lcomp)
                  error('index cannot exceed the number of selected components')
              end
              comp1 = comp(indc);
              lcomp1 = length(comp1);
              Cc1 = Cc(:,indc);
              
              if nargs == 2
                 for i = 1:lcomp1
                     figure 
                     tsplot(Cc1(:,i),datei,comp1{i})
                     pause
                 end
              end
             
          elseif (iscellstr(varargin{2})) && (nargs == 2) 
              
             vnames = varargin{2};
             
             if ~isequal(length(vnames),lcomp)
                 error('number of components must agree with the number of names')
             end
             for i = 1:lcomp
                 figure 
                 tsplot(Cc(:,i),datei,vnames{i})
                 pause
             end
          else 
              error('fifth argument must be a scalar or a cell with scalars if there are more than five inputs to dispcomp')
          end
   
       
         if nargs == 3
             
             if (~iscellstr(varargin{3})) && (~ischar(varargin{3}))
                error('sixth argument must be a string or string array') 
             end
            
             vnames = varargin{3};
             vnames = {vnames};
             
             if ~isequal(length(vnames),lcomp1)
                 error('number of components to be plotted must agree with the number of names')
             end
             
             for i = 1:lcomp1
                 figure 
                 tsplot(Cc1(:,i),datei,vnames{i})
                 pause
             end
  
          
         elseif nargs > 3
             error('maximum number of arguments is six')
         end       
   end    
end
   
end