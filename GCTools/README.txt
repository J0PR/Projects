adjMatrix=causalMatrix(data(nVars,nObs),Measure,Model,type,Fs,fArray,MINP,MAXP,scrit,pval,nullPopMethod,numSurrogates,nullStruct);
or time-varying
adjMatrix=TV_CausalMatrix(data(nVars,nObs),Measure,Model,type,Fs,fArray,MINP,MAXP,scrit,pval,nullPopMethod,numSurrogates,nullStruct);

Before using any function please add the folder functions (and its subfolders) to your workspace.

For causalMatrix()

data is in nº variables vs nº observations

Measure:
	'GCGeweke'
	'DTF'
	'PDC'
	'DC'
	'PartitionsGCGeweke'
	'GCI'
	'newGCTime'
	'TE'
	'PSI'

Model
	'LS'	stands for Least squares
	'MVAR'	stands for multivariate vector AR
	'CP'	stands for correlation purged Least squares
	'NP'	stands for non-parametric
	'knn'   stands for k-nearest neighbours. To be used only with TE measure!

type
	'pwise'	stands for pairwise analysis
	'cond'	stands for conditional analysis

Fs	sampling frequency in Hz

fArray	frequncy range array in Hz. 
For example 
	to analyse frequency between 5 and 30 Hz: fArray=[5 30]; returns a adjMatrix of dim nVars*nVars
	to get several (in the exmple 3) frequency bands (5->10Hz, 10->20Hz,50->100Hz): fArray=[5,10;10,20;50,100]; returns a adjMatrix of dim nVars*nVars*3
	to get all frequncy bins: fArray = 'all' returns a adjMatrix of dim nVars*nVars*length(all frequencies from 0 to nyquist freq)/df
	to get an average off all frequencies (no frequency discrimination): fArray='mean';

MINP,MAXP
	min and maximum model orderl. Search will be made with BIC or AIC according to the scrit variable
	if MINP==MAXP the orther is fixed with their value.
	
scrit
	'bic' Bayesian information criteria
	'aic' Akaike information criteria

pval
	pvalue for significance testing

nullPopMethod
	method to generate surrogates for null hypothesis testing. Use 'phaseran' or 'AAFT'.

numSurrogates
	Number of surrogates. 100 is a good value.. 0 avoids statistical significance

nullStruct
	keep it empty [] so significance testing is done. 
	


EXAMPLE:
GCCond=causalMatrix(rand(3,100),'GCI','LS','cond',250,[5 30],1,30,'bic',0.01,'phaseran',100,[]);

%%%_________________________________________________________________________________________________________%%%

For TV_CausalMatrix()

data is in nº variables vs nº observations

Measure:
	'GCGeweke'
	'DTF'
	'PDC'
	'DC'
	'PartitionsGCGeweke'
	'GCI'
	'newGCTime'

Model
	'TV_MVAR'	time-varying stands for multivariate vector AR
	'TV_VARMA'	time-varying stands for multivariate vector ARMA

type
	'pwise'	stands for pairwise analysis
	'cond'	stands for conditional analysis

Fs	sampling frequency in Hz

fArray	frequncy range array in Hz. 
For example 
	

MINP,MAXP
	min and maximum model orderl. Search will be made with BIC or AIC according to the scrit variable
	if MINP==MAXP the orther is fixed with their value.
	
scrit
	'bic' Bayesian information criteria
	'aic' Akaike information criteria

pval
	pvalue for significance testing

nullPopMethod
	Not working yet in TV_CausMat!!!

numSurrogates
	Not working yet in TV_CausMat!!!

nullStruct
	keep it empty [] so significance testing is done. 
	


EXAMPLE:
GCCond=causalMatrix(rand(3,100),'GCI','LS','cond',250,[5 30],1,30,'bic',0.01,'phaseran',100,[]);