//####################################################################
//Mutual information estimator
//####################################################################
//05 July 2013
//
//####################################################################

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "mex.h"

int neiE(double **x, int i, int comp1, int comp2, int dim, int bs, double epsgrid, double eps, int **box, int *lis) {
  int ix,iy,ix1,iy1,ix2,jj,step,d;
  int el,nx,ib=bs-1;
  double dd,dy;
  double *xx;

  xx=(double *)mxCalloc(dim,sizeof(double));
  for (d=0;d<dim;d++) xx[d]=x[d][i];

  ix=(int)(xx[comp1]/epsgrid)&ib; iy=(int)(xx[comp2]/epsgrid)&ib;
  jj=0; nx=0;
  while (eps>epsgrid*(jj-1)) {
    step = (jj) ? 2*jj : 1;
    for (ix1=ix-jj;ix1<=ix+jj;ix1++) {
      ix2=ix1&ib;
      for (iy1=iy-jj;iy1<=iy+jj;iy1+=step) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	  if (dd<=eps) nx++;
	  el=lis[el];
	}
      }
    }
    for (ix1=ix-jj;ix1<=ix+jj;ix1+=step) {
      ix2=ix1&ib;
      for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	  if (dd<=eps) nx++;
	  el=lis[el];
	}
      }
    }
    jj++;
    if (jj==(bs/2)) break;
  }
  if ( jj==(bs/2) ) { //half of the layer
    for (ix1=ix-jj;ix1<ix+jj;ix1++) {
      ix2=ix1&ib; iy1=iy-jj;
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	if (dd<=eps) nx++;
	el=lis[el];
      }
    }
    ix1=ix-jj; ix2=ix1&ib;
    for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	if (dd<=eps) nx++;
	el=lis[el];
      }
    }
  }
  mxFree(xx);
  return nx-1;
}

int neiE1(double *x, int i, double scal, int bs, double eps,  int *box, int *lis, int *mxi) {
  double dd;
  int mi,mp,mm,nx=0;
  double xc=x[i];

  mp=(int)((xc+eps)*scal); if(mp>bs) mp=bs;
  mm=int((xc-eps)*scal); if(mm<0)  mm=0;
  mi=box[mp];
  while(mi>=0) {
    dd=x[mi]-xc;if(fabs(dd)<=eps) nx++;
    mi=lis[mi];
  }
  if(mm>=mp) return nx-1;
  mi=box[mm];
  while(mi>=0) {
    dd=xc-x[mi];if(fabs(dd)<=eps) nx++;
    mi=lis[mi];
  }
  nx+=mxi[mp-1]-mxi[mm];
  return nx-1;
}

void neiK(double **x, int dim, int comp1, int comp2, int i, 
	  int bs, double epsgrid, int K, int **box, int *lis,
	  int *nn) {
  double *dn,*xx;
  int k,ix,iy,ix1,iy1,ix2,jj,step,ib=bs-1;
  int el;
  double dd,dy;
  int d;

  dn=(double*)mxCalloc(K+1,sizeof(double));
  xx=(double*)mxCalloc(dim,sizeof(double));
  for(k=0;k<dim;k++) xx[k]=x[k][i];
  dn[0]=0;
  for(k=1;k<=K;k++) dn[k]=1e30;
  ix=(int)(xx[comp1]/epsgrid)&ib; iy=(int)(xx[comp2]/epsgrid)&ib;
  jj=0; 
  while (dn[K]>epsgrid*(jj-1)) {
    step = (jj) ? 2*jj : 1;
    for (ix1=ix-jj;ix1<=ix+jj;ix1++) {
      ix2=ix1&ib;
      for (iy1=iy-jj;iy1<=iy+jj;iy1+=step) {
	el=box[ix2][iy1&ib];
	while (el != -1) { if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	} el=lis[el];
	}
      }
    }
    for (ix1=ix-jj;ix1<=ix+jj;ix1+=step) {
      ix2=ix1&ib;
      for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
	el=box[ix2][iy1&ib];
	while (el != -1) { if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if ( (dy=fabs(xx[1]-x[1][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  } 
	} el=lis[el];
	}
      }
    }
    jj++;
    if (jj==bs/2) break;
  }

  if ( jj==(bs/2) ) { //half of the layer
    for (ix1=ix-jj;ix1<ix+jj;ix1++) {
      ix2=ix1&ib; iy1=iy-jj;
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	}
	el=lis[el];
      }
    }
    ix1=ix-jj; ix2=ix1&ib;
    for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	} 
	el=lis[el];
      }
    }
  }
  mxFree(dn);mxFree(xx);
}
void make_box1(double *x, int N, double scal, int bs, 
	       int *box, int *lis, int *mxi) {
  int i, ix;
  for (i=0;i<=bs;i++) box[i]=-1;
  for (i=0;i<=bs;i++) mxi[i]=0;
  for(i=0;i<N;i++) {ix=(int)(x[i]*scal);lis[i]=box[ix]; box[ix]=i; mxi[ix]++;}
  for (i=1;i<=bs;i++) mxi[i]+=mxi[i-1];
}

void make_box2ind(double **x, int dim, int N, int comp1, int comp2, int bs, int inveps, 
		  int *ind,  int **box, int *lis) {
  int d,i,ix,iy,ixy;
  int ib=bs-1;
  double **xx;

  xx=(double **)mxCalloc(dim, sizeof(double*)); 
  for (d=0;d<dim;d++) {
    xx[d]=(double *)mxCalloc(N, sizeof(double)); 
    memcpy(xx[d],x[d],N*sizeof(double));
  }
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) box[ix][iy] = -1;
  
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  i=-1;
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) {
    ixy=box[ix][iy];
    while(ixy>=0) {
      i++;
      for (d=0;d<dim;d++) x[d][i]=xx[d][ixy];
      ind[ixy]=i;ixy=lis[ixy];
    }
    box[ix][iy]=-1;
  }
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  for (d=0;d<dim;d++) mxFree(xx[d]); mxFree(xx);
}
void mir_xnyn(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mir) {
  int i,k,nx2,ny2;
  double *xc,dy,dx;
  double epsx,epsy;
  double dxy2;
  double **xx,**yy;;
  double scalx, scaly;
  int *nn;
  int d;

  int BOX,BOX1;
  int **box,**boxy,**boxx,*lis; // two dimensional boxes
  int *lisy,*lisx; // lists for two dimensions
  int *boxx1, *boxy1; // onedimensional boxes
  int *lisy1,*lisx1; // lists for one dimensions
  int *mxi, *myi; //accumulative lists of points in oned boxes
  int *ind,*indx,*indy; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)mxCalloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-1/(double(i));   // 

  nn=(int*)mxCalloc(K+1,sizeof(int));

  xc=(double*)mxCalloc(dimx+dimy,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  if (dimx>1) {
    xx=(double**)mxCalloc(dimx,sizeof(double*));
    for (d=0;d<dimx;d++) xx[d]=(double*)mxCalloc(N,sizeof(double));
    boxx=(int**)mxCalloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxx[i]=(int*)mxCalloc(BOX,sizeof(int));
    lisx=(int*)mxCalloc(N,sizeof(int));
  } else { boxx1=(int*)mxCalloc(BOX1+1,sizeof(int)); 
  lisx1=(int*)mxCalloc(N,sizeof(int)); mxi=(int*)mxCalloc(BOX1+1,sizeof(int)); }
  if (dimy>1) {
    yy=(double**)mxCalloc(dimy,sizeof(double*));
    for (d=0;d<dimy;d++) yy[d]=(double*)mxCalloc(N,sizeof(double));
    boxy=(int**)mxCalloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxy[i]=(int*)mxCalloc(BOX,sizeof(int));
    lisy=(int*)mxCalloc(N,sizeof(int));
  } else { boxy1=(int*)mxCalloc(BOX1+1,sizeof(int)); 
  lisy1=(int*)mxCalloc(N,sizeof(int)); myi=(int*)mxCalloc(BOX1+1,sizeof(int)); }
 
  box=(int**)mxCalloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)mxCalloc(BOX,sizeof(int));
  lis=(int*)mxCalloc(N,sizeof(int));

  ind=(int*)mxCalloc(N,sizeof(int));
  indx=(int*)mxCalloc(N,sizeof(int));
  indy=(int*)mxCalloc(N,sizeof(int));

  //save x if it would be reordered
  if (dimx>1) for (d=0;d<dimx;d++) memcpy(xx[d],x[d],N*sizeof(double)); 
  if (dimy>1) for (d=0;d<dimy;d++) memcpy(yy[d],x[d+dimx],N*sizeof(double));

  make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps,ind,box,lis); 
  //for searching neighbours in product space

  if (dimx==1) {scalx=scal[0]; make_box1(x[0],N,scalx,BOX1,boxx1,lisx1,mxi);}
  else make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps,indx,boxx,lisx); 
  if (dimy==1) {scaly=scal[dimx]; make_box1(x[dimx],N,scaly,BOX1,boxy1,lisy1,myi); }
  else make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps,indy,boxy,lisy); 

  dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dimx+dimy;d++) xc[d]=x[d][ind[i]];

    neiK(x,dimx+dimy,0,dimx,ind[i],BOX,epsilon,K,box,lis,nn);
    epsx=0; for (d=0;d<dimx;d++) 
      for(k=1;k<=K;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx ) epsx=dx;
    epsy=0; for (d=dimx;d<dimx+dimy;d++) 
      for(k=1;k<=K;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy ) epsy=dy;
    if (dimx>1) nx2=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,epsx,boxx,lisx);
    else  nx2=neiE1(x[0],ind[i],scalx,BOX1,epsx,boxx1,lisx1,mxi);
    if (dimy>1) ny2=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,epsy,boxy,lisy);
    else ny2=neiE1(x[dimx],ind[i],scaly,BOX1,epsy,boxy1,lisy1,myi);

    dxy2+=psi[nx2]+psi[ny2];
  }
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;

  mxFree(xc);mxFree(nn);mxFree(lis);
  for (i=0;i<BOX;i++) mxFree(box[i]); mxFree(box);
  mxFree(ind);mxFree(indx);mxFree(indy);
  if (dimx==1) {mxFree(mxi);mxFree(boxx1);mxFree(lisx1);} 
  else { for (i=0;i<BOX;i++) mxFree(boxx[i]); mxFree(boxx); mxFree(lisx); for (d=0;d<dimx;d++) mxFree(xx[d]); mxFree(xx); }
  if (dimy==1) {mxFree(myi);mxFree(boxy1);mxFree(lisy1);} 
  else { for (i=0;i<BOX;i++) mxFree(boxy[i]); mxFree(boxy); mxFree(lisy); for (d=0;d<dimy;d++) mxFree(yy[d]); mxFree(yy); }
  mxFree(phi);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      /* Macros for the ouput and input arguments */
#define B_OUT plhs[0]
#define X_IN prhs[0]
#define NdX_IN prhs[1]
#define NdY_IN prhs[2]
#define N_IN prhs[3]
#define K_IN prhs[4]
    

  int i;
  double **x;
  double *scal;
  double *min;
  double *max;
  double *psi;
  int K,N;
  int d;
  double mir;
  int dimx,dimy;


  int BOX1;

  double s,me;
  double addnoise=-1;


double *B, *A;
int M, NN, m, n;

dimx = (int)mxGetScalar(NdX_IN); /* Get dimx */
dimy = (int)mxGetScalar(NdY_IN); /* Get dimy */
N = (int)mxGetScalar(N_IN); /* Get N */
K = (int)mxGetScalar(K_IN); /* Get K */

M = (int)mxGetM(X_IN); /* Get the dimensions of X */
NN = (int)mxGetN(X_IN);
A = mxGetPr(X_IN); /* Get the pointer to the data of X_IN */
B_OUT = mxCreateDoubleMatrix(1, 1, mxREAL); /* Create the output matrix */
B = mxGetPr(B_OUT); /* Get the pointer to the data of B */

  x=(double**)mxCalloc(dimx+dimy,sizeof(double*));
  for (d=0;d<dimx+dimy;d++) x[d]=(double*)mxCalloc(N,sizeof(double));
  scal=(double*)mxCalloc(dimx+dimy,sizeof(double));
  min=(double*)mxCalloc(dimx+dimy,sizeof(double));
  max=(double*)mxCalloc(dimx+dimy,sizeof(double));
  for (d=0;d<dimx+dimy;d++) {min[d]=DBL_MAX/2;max[d]=-DBL_MAX/2;}
  
  for(n = 0; n < NN; n++) /* reshape input matrix */
      for(m = 0; m < M; m++) x[m][n] = A[m + M*n];
  
  // add noise
  if (addnoise) {
    srand((dimx+dimy)*N*K*(int)(x[(dimx+dimy)/2][N/10]));
    if (addnoise==-1) for (d=0;d<dimx+dimy;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*1e-8;
    else for (d=0;d<dimx+dimy;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*addnoise;
  }

  //normalization
  for (d=0;d<dimx+dimy;d++) {
    me=s=0; for (i=0;i<N;i++) me+=x[d][i];
    me/=N;  for (i=0;i<N;i++) s+=(x[d][i]-me)*(x[d][i]-me);
    s/=(N-1);s=sqrt(s);
    if (s==0) {;}
    for (i=0;i<N;i++) {
      x[d][i] = (x[d][i]-me)/s;
      if (x[d][i]<min[d]) min[d]=x[d][i]; 
      if (x[d][i]>max[d]) max[d]=x[d][i];
    }
    for (i=0;i<N;i++) x[d][i]=x[d][i]-min[d];
  }

  psi=(double*)mxCalloc(N+1,sizeof(double)); 
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i;
  BOX1=N-5;
  for (d=0;d<dimx+dimy;d++) scal[d]=BOX1/(max[d]-min[d]); 

  mir_xnyn(x,dimx,dimy,N,K,psi,scal,&mir);
  *B=mir;


  for (d=0;d<dimx+dimy;d++) mxFree(x[d]); mxFree(x);
  mxFree(scal);
  mxFree(min);mxFree(max);
  mxFree(psi);
  return;
}

