#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include "memory.h"
#include <math.h> 
#include <R.h>
#include <Rmath.h>
#include <S.h>

#define Integer long
#define Real double



Real geoRmatern(Real uphi, Real kappa)
{   
  
  /* 
     WARNING: THIS FUNCTION IS COPIED IN geoRglmm
     NOTIFY OLE ABOUT ANY CHANGE 
  */
  
  Real ans,cte;
  
  if (uphi==0) return 1;
  else{
    if (kappa==0.5) 
      ans = exp(-uphi);
    else {
      cte = pow(2, (-(kappa-1)))/gammafn(kappa); 
      ans = cte * R_pow(uphi, kappa) * bessel_k(uphi,kappa,1); 
    }
  }
  /* Rprintf("   ans=%d ", ans); */
  return ans; 
}

Real corrfctvalue(Real phi, Real kappa, Real h, Integer cornr)
{
  
  /* 
     WARNING: THIS FUNCTION IS COPIED IN geoRglmm
     NOTIFY OLE ABOUT ANY CHANGE 
     
     Correlation functions implemented and their numbers
     
     1: PURE NUGGET
     2: EXPONENTIAL
     3: SPHERICAL
     4: GAUSSIAN
     5: WAVE (hole effect)
     6: CUBIC
     7: POWER
     8: POWERED EXPONENTIAL
     9: CAUCHY
     10: GNEITING
     11: CIRCULAR
     12: MATERN
     13: GNEITING-MATERN (NOT YET IMPLEMENTED)
     
     WARNING: codes above must be the same as in the geoR/geoS function
     "cor.number"
  */
  
  Real hphi, hphi2, hphi4;
  if(h==0) return 1;
  else{  
    hphi  = h/phi ;
    switch(cornr){
    case 1: /* pure nugget */
      return 0 ;
      break;
    case 2: /* exponential */
      return exp(-hphi) ;
      break;
    case 3: /* spherical */
      if (h < phi) 
	return 1 - (1.5 * hphi) + (0.5 * hphi*hphi*hphi) ;
      else
	return 0 ;
      break;
    case 4: /* Gaussian */
      return exp(-(hphi * hphi)) ;
      break;
    case 5: /* wave (hole effect) */
      return hphi*sin(hphi) ;
      break;
    case 6: /* cubic */
      if (h < phi){
	hphi2 = hphi * hphi ;
	hphi4 = hphi2 * hphi2 ;
	return 1 - ((7 * hphi2) - (8.75 * hphi2 * hphi) + (3.5 * hphi4 * hphi) - (0.75 * hphi4 * hphi2 * hphi)) ;
      }
      else
	return 0 ;
      break;
    case 7: /* power */
      return exp(phi * log(h)) ;
      break;
    case 8: /* powered.exponential */
      return exp(-1 *  R_pow(hphi, kappa))  ;
      break;
    case 9:  /* cauchy */
      return R_pow((1 + (hphi * hphi)), (-kappa)) ;
      break;
    case 10:  /* gneiting */
      hphi4 = 1 - hphi;
      if (hphi4 > 0) hphi4 = pow(hphi4, 8);
      else hphi4 = 0 ;
      hphi2 = hphi * hphi ;
      return (1 + 8 * hphi + 25 * hphi2 + 32 * (hphi2 * hphi)) * hphi4 ;
      break;
    case 11: /* circular */
      if(h < phi){
	return  1 - (M_1_PI * (2 * (hphi * sqrt(1 - hphi * hphi)
				    + asin(hphi)))) ;
      }
      else
	return 0 ;
      break;
    case 12: /* matern */
      return geoRmatern(hphi, kappa);
      break;
      /* case 13: gneiting-matern NOT YET IMPLEMENTED 
	 res[ind] =  ;
	 break; */
    default: 
      return -1;
      break;
    }
  }
}

void veccorrval(Real *phi, Real *kappa, Real *h, Integer *n, 
		Integer *cornr, Real *res)  
{
  Integer register j ;
  
  for (j=0; j<(*n); j++)
    res[j] = corrfctvalue(*phi, *kappa, h[j], *cornr)  ;
  
  return ;
}

void diag_quadraticform_XAX(Real *lower, Real *diag, Real *xvec, 
			    Integer *nx, Integer *n, Real *res)
     /*
       This function computes quadratic forms of the type x'Ax
       where A is a symmetric matrix.
       If X is a matrix (passed to C in a vector form) the quadratic
       form is computed for each column of X
       
       lower  : lower triangle of the A matrix
       diag   : diagonal of the A matrix
       xvec   : the vector(s) for the quad. form.
       nx     : number of vectors of type x 
       (= number of quadratic forms to be computed)
       n      : size of each vector
       res    : store the result(s)
       
       author : Paulo J.Ribeiro Jr , 03/11/00
     */
     
{
  Integer pos;
  Integer register i,j,k;
  Real xij, xii;
  
  for(k=0; k<*nx; k++){
    
    pos = 0;
    xij = 0.0;
    for (j=0; j<(*n-1); j++) {
      for (i=j+1; i<*n; i++) {
	xij += (xvec[(k*(*n) + j)] * xvec[(k*(*n) + i)] * lower[pos++]) ;
      }
    }
    
    pos = 0;
    xii = 0.0;
    for (i=0; i<*n; i++){
      xii += (xvec[(k*(*n) + i)] * xvec[(k*(*n) + i)] * diag[i]) ;
    }
    res[k] = xii + 2*xij ;    
  }
  
}

void bilinearform_XAY(Real *lower, Real *diag, Real *xvec, 
		       Real *yvec, Integer *nx, Integer *ny, 
		       Integer *n, Real *res)
     /*
       This function computes retangular forms of the type x'Ay
       where A is a symmetric matrix.
       X and Y can be matrices in whic cases quadratic formas are computed
       for each row of X combined with each column of Y
       
       lower  : lower triangle of the A matrix
       diag   : diagonal of the A matrix
       xvec   : the vector(s) for the quad. form.
       yvec   : the vector y in the quad. form.
       nx     : number of vectors of type x
       (= number of quadratic forms to be computed)
       ny     : number of vectors in yvec
       n      : dimension of the matix A 
       (also size of each vector yvec and xvec)
       res    : store the result(s)
       
       author: Paulo J.Ribeiro Jr , 03/11/00
     */
     
{
  Integer pos;  
  Integer register i,j,k,l;
  Real xyij, xyji, xyii;
  
  for (l=0; l<*ny; l++){
    
    for(k=0; k<*nx; k++){
      pos = 0;
      xyij = 0.0;
      xyji = 0.0;
      for (j=0; j<(*n-1); j++) {
	for (i=j+1; i<*n; i++) {
	  xyji += (xvec[(k*(*n) + j)] * lower[pos] * yvec[(l*(*n) + i)]) ;
	  xyij += (xvec[(k*(*n) + i)] * lower[pos] * yvec[(l*(*n) + j)]) ;
	  pos++ ;
	}
      }
      xyii = 0.0;
      for (i=0; i<*n; i++){
	xyii += (xvec[(k*(*n) + i)] * diag[i] * yvec[(l*(*n) + i)]) ;
      }
      res[(l * (*nx) + k)] = xyii + xyij + xyji ;
    }
  }
}

void loccoords(Real *xloc, Real *yloc, Real *xcoord, Real *ycoord, 
	       Integer *nl, Integer *nc, Real *res) 
     /* This function computes the distance between each data location
	to each of the prediction location
	
        xloc, yloc     : xy coordinates of the prediction locations
        xcoord, ycoord : xy coordinates of the data points
	nl, nc         : number of prediction locations and data locations
	res            : stores the results to be returned, 
	a vector with distances
     */   
     
{ 
  Integer register i,j, ind;
  Real register dx,dy;
  
  ind = 0;
  for (j=0; j<*nl; j++) {  
    for (i=0; i<*nc; i++) {
      dx = (xloc[j] - xcoord[i]) ;
      dy = (yloc[j] - ycoord[i]) ;
      res[ind++] = pythag(dx,dy) ;
    }
  }
  
}

void tgangle(Real *xloc, Real *yloc, Integer *nl, Real *res) 
     /* 
	This function computes the tangent of the (azimuth) 
        angle between pairs of locations
	
        xloc, yloc     : xy coordinates of the locations
	nl,            : number of locations
	res            : stores the results to be returned, 
	                 a vector with tangent of the angles
     */   
     
{ 
  Integer register i,j, ind;
  Real register dx,dy;
  
  ind = 0;
  for (j=0; j<*nl; j++) {  
    for (i=j+1; i<*nl; i++) {
      dx = (xloc[i] - xloc[j]) ;
      dy = (yloc[i] - yloc[j]) ;
      res[ind] = dy/dx ;
      ind++ ;
    }
  }
  
}

void distdiag(Real *xloc, Real *yloc, Integer *nl, Real *res) 
     /* This function computes the distances between locations
	including the diagonal term
	
        xloc, yloc     : xy coordinates of the locations
	nl,            : number of locations
	res            : stores the results to be returned, 
	                 a vector with distances
     */   
     
{ 
  Integer register i,j, ind;
  Real register dx,dy;
  
  ind = 0;
  for (j=0; j<*nl; j++) {  
    for (i=j; i<*nl; i++) {
      if(i==j)
	res[ind] = 0.0 ;
      else{
	dx = (xloc[j] - xloc[i]) ;
	dy = (yloc[j] - yloc[i]) ;
	res[ind] = pythag(dx, dy) ;
      }
      ind++ ;
    }
  }
  
}

void binit(Integer *n, Real *xc, Real *yc, Real *sim, 
	   Integer *nbins, Real *lims, Integer *robust, 
	   Real *maxdist, Integer *cbin, Real *vbin,
	   Integer *sdcalc, Real *sdbin)
{
  
  Integer register i, j, ind=0;
  Real register v=0.0;
  Real dist=0.0, dx=0.0, dy=0.0;
  
  for (j=0; j < *n; j++)
    { 
      for (i=j+1; i<*n; i++) 
	{
	  dx = xc[i] - xc[j];
	  dy = yc[i] - yc[j];
	  dist = pythag(dx, dy);
	  
	  if(dist <= *maxdist)
	    {
	      v = sim[i] - sim[j];
	      if (*robust) v = sqrt(sqrt(v*v));
	      else v = (v*v)/2.0;
	      ind = 0;
	      while (dist >= lims[ind] && ind <= *nbins ) ind++ ;
	      if (dist < lims[ind])
		{
		  vbin[(ind-1)]+= v; 
		  cbin[(ind-1)]++;
		  if(*sdcalc) sdbin[(ind-1)] += v*v;
		}
	    }
	}
    }
  
  for (j=0; j < *nbins; j++) 
    {
      if (cbin[j]){
	if(*sdcalc)
	  { 
	    sdbin[j] = sqrt((sdbin[j] - ((vbin[j] * vbin[j])/cbin[j]))/(cbin[j] - 1));
	  }
	vbin[j] = vbin[j]/cbin[j];
	if (*robust) {
	  vbin[j] = vbin[j] * vbin[j];
	  vbin[j] = (vbin[j] * vbin[j])/(0.914 + (0.988/cbin[j]));
	}
      }
    }
  
}

     
void lower_DIAGminusXAX(Real *lower, Real *diag, Real *xvec, 
	      Integer *nxcol, Integer *n, Real *Dval, Real *res)
     /*
       This function computes the lower triangle and diagonal
       of forms (D - X'AX), where D is a diagonal matrix, 
       A is a symmetric matrix and X a retangular matrix.
              
       lower  : lower triangle of the A matrix
       diag   : diagonal of the A matrix
       xvec   : matrix X in a vector form
       nxcol  : number of coluns of matrix X
       n      : dimension of the matix A 
       Dval   : element of the diagonal matrix
       res    : store the result
       
       author : Paulo J.Ribeiro Jr , 01/12/00
     */
     
{
  Integer pos;  
  Integer register i,j,k,l;
  Real xyij, xyji, xyii;
  
  for (l=0; l<*nxcol; l++){
    
    for(k=l; k<*nxcol; k++){
      pos = 0;
      xyij = 0.0;
      xyji = 0.0;
      for (j=0; j<(*n-1); j++) {
	for (i=j+1; i<*n; i++) {
	  xyji += (xvec[(k*(*n) + j)] * lower[pos] * xvec[(l*(*n) + i)]) ;
	  xyij += (xvec[(k*(*n) + i)] * lower[pos] * xvec[(l*(*n) + j)]) ;
	  pos++ ;
	}
      }
      xyii = 0.0;
      for (i=0; i<*n; i++){
	xyii += (xvec[(k*(*n) + i)] * diag[i] * xvec[(l*(*n) + i)]) ;
      }
      if (k > l)
	res[(*nxcol * l - (l * (l+1)/2) + k)] = -1 * (xyii + xyij + xyji) ;
      else
	res[(*nxcol * l - (l * (l+1)/2) + k)] = (*Dval) - (xyii + xyij + xyji) ;
    }
  }
}

void lower_R0minusXAXplusBvar(Real *lower, Real *diag, Real *xvec, 
			      Integer nxcol, Integer n, Real *Dval, 
			      Real *Blower, Real *Bdiag, Real *bvec, 
			      Integer Bsize, Real *ss, Real *res)
     /*
       This function computes the lower triangle and diagonal
       of forms (R0 - X'AX + b'Bb), where R0 is the cov.  matrix
       between prediction points; 
       A is a symmetric matrix and X a retangular matrix corresponding
       to the reduction in variance due to the data;
       B is symmetric and b is a retaulkar matrix correspond in the 
       increase in the variance due to the unknown mean
              
       lower  : lower triangle of the A matrix
       diag   : diagonal of the A matrix
       xvec   : matrix X in a vector form
       nxcol  : number of coluns of matrix X
       n      : dimension of the matrix A 
       Dval   : element of the diagonal matrix
       Blower : lower triangle  of the Var(beta) matrix 
       Bdiag  : diagonal  of the Var(beta) matrix 
       bvec   : t(b) matrix in vector form
       Bsize  : dimension of mean vector beta and its variance matrix
       ss     : 
       res    : input lower triangle of R0 and store the result

       
       author : Paulo J.Ribeiro Jr , 04/12/00
     */
     
{
  Integer pos, bpos, indpos=0;  
  Integer register i,j,k,l;
  Real xyij, xyji, xyii;
  Real bxyij, bxyji, bxyii;
  
  for (l=0; l<nxcol; l++){
    
    for(k=l; k<nxcol; k++){
      
      /*      
	      Computing lower triangle (including diagonal) of XAX
      */
      
      pos = 0;
      xyij = 0.0;
      xyji = 0.0;
      for (j=0; j<(n-1); j++) {
	for (i=j+1; i<n; i++) {
	  xyji += (xvec[(k * n + j)] * lower[pos] * xvec[(l * n + i)]) ;
	  xyij += (xvec[(k * n + i)] * lower[pos] * xvec[(l * n + j)]) ;
	  pos++ ;
	}
      }
      xyii = 0.0;
      for (i=0; i<n; i++){
	xyii += (xvec[(k * n + i)] * diag[i] * xvec[(l * n + i)]) ;
      }
      
      /*
	Computing lower triangle (including diagonal) of bBb
      */

      bpos = 0;
      bxyij = 0.0;
      bxyji = 0.0;
      bxyii = 0.0;
      if (Bsize == 1){
	bxyii = bvec[l] * bvec[k] * (*Bdiag) ;
      }      
      else{
	for (j=0; j<(Bsize-1); j++) {
	  for (i=j+1; i<Bsize; i++) {
	    bxyji += (bvec[(k*(Bsize) + j)] * Blower[bpos] * bvec[(l*(Bsize) + i)]) ;
	    bxyij += (bvec[(k*(Bsize) + i)] * Blower[bpos] * bvec[(l*(Bsize) + j)]) ;
	    bpos++ ;
	  }
	}
	for (i=0; i<Bsize; i++){
	  bxyii += (bvec[(k*(Bsize) + i)] * Bdiag[i] * bvec[(l*(Bsize) + i)]) ;
	}
      }
      
      /*      
	      Computing lower triangle (including diagonal) of (R0 - XAX + bBb)
	      indpos = (nxcol * l - (l * (l+1)/2) + k) ;
	      Rprintf("\n indpos=%d ", indpos);
      */
      
      if (k > l){
	res[indpos] += (-1 * (xyii + xyij + xyji) + (bxyii + bxyij + bxyji)) ;
	res[indpos] *= *ss ;
      }
      else{
	res[indpos] *= (*Dval) ;
	res[indpos] += (-1 * (xyii + xyij + xyji)  + (bxyii + bxyij + bxyji)) ;
	res[indpos] *= *ss ;
      }
      indpos++ ;
    }
  }
}


void chol(Real *inmatrix, Integer N)
{
  Integer register Drow, Dcol, Dcol2;
  Real register sum;
  Real register *Pcol, *anothercol;
  /* 
     returns L where L L'=inmatrix
     NR function choldc, sec. 2.9
     i=dcol, j=drow, k=dcol2
  */
  for(Dcol=0;Dcol<N;Dcol++) {
    Pcol=inmatrix + Dcol * N - ((Dcol * (Dcol+1))/2);

    for(Drow=Dcol;Drow<N;Drow++){
      for(sum=Pcol[Drow],Dcol2=Dcol-1;Dcol2>=0;Dcol2--) {
	anothercol = inmatrix + N * Dcol2 - ((Dcol2 * (Dcol2+1))/2);
	/*
	sum -= inmatrix[matref(Drow,Dcol2,N)]*inmatrix[matref(Dcol,Dcol2,N)];
	*/	
	sum -= anothercol[Drow]*anothercol[Dcol];
      }
      if (Drow == Dcol) {
	if (sum<=0.0) {
	  error("%s%ld%s%e", "chol: matrix not pos def, diag[" , Drow , "]= " , sum);
	  return;
	}
	Pcol[Drow]=sqrt(sum);
      } else Pcol[Drow]=sum/Pcol[Dcol];
    }
  }
  /* 
     at this point the diagonal and lower triangle of inmatrix
     contain the cholesky decomp 
  */
  return;
}

void mvnorm(Real *means, Real *Q, Real *nscores, Integer N, 
	    Integer Nsims, Real *Vsqglchi) 
{  
  /*  
      returns means + chol(Q) %*% nscores 
  */
  Integer register Drow, Dcol, Dsim, i ;
  Real *Vsim = (Real*) malloc(sizeof(Real)*N) ;
  
  chol(Q, N);  
  
  /*
    multiplyLower(X, Q, nscores, *N);
  */

  /* this was wrong before
    for(Dsim=0;Dsim<Nsims; Dsim++){    
    for(Drow=0;Drow<N;++Drow) {
    Vsim[Drow] = means[Drow] ;
    for(Dcol=0;Dcol<=Drow;++Dcol) {
    Vsim[Drow] += Q[ N * Dcol - ((Dcol * (Dcol+1))/2) + Drow] * nscores[((N * Dsim) + Dcol)];
    }
    }
    for(i=0;i<N;i++){
    nscores[(((N) * Dsim) + i)] = Vsim[i] * Vsqglchi[Dsim];
    }
    }
  */
  
  for(Dsim=0;Dsim<Nsims; Dsim++){    
    for(Drow=0;Drow<N;++Drow) {
      Vsim[Drow] = 0.0 ;
      for(Dcol=0;Dcol<=Drow;++Dcol) {
	Vsim[Drow] += Q[ N * Dcol - ((Dcol * (Dcol+1))/2) + Drow] * nscores[((N * Dsim) + Dcol)];
      }
    }
    for(i=0;i<N;i++){
      nscores[(((N) * Dsim) + i)] = means[i] + (Vsim[i] * Vsqglchi[Dsim]);
    }
  }
  free(Vsim);
  
}

void multmvnorm(Real *means, Real *Q, Real *nscores, Integer N, 
		Integer Nsims, Real *Vsqglchi) 
{
  /* 
     returns means + chol(Q) %*% nscores, 
     where means is a matrix (i.e. means can be different for each simulation) 
  */ 
  
  Integer Drow, Dcol, Dsim ;
  Real *Vsim = (Real*) malloc(sizeof(Real)*N) ;
  chol(Q, N);  
  for(Dsim=0;Dsim<Nsims; Dsim++){    
    for(Drow=0;Drow<N;++Drow) {
      Vsim[Drow] = 0.0 ;
      for(Dcol=0;Dcol<=Drow;++Dcol) {
	Vsim[Drow] += Q[ N * Dcol - ((Dcol * (Dcol+1))/2) + Drow] * nscores[((N * Dsim) + Dcol)];
      }
    }
    for(Drow=0;Drow<N;++Drow) {
      nscores[(((N) * Dsim) + Drow)] = means[(N * Dsim)+Drow] + (Vsim[Drow] * Vsqglchi[Dsim]);
    }
  }
}

void multiplyLower(Real *X, Real *A, Real *B, Integer *N) 
{  
  Integer register Drow, Dcol;
  
  for(Drow=0;Drow<*N;++Drow) {
    X[Drow]=0;
    for(Dcol=0;Dcol<=Drow;++Dcol) {
      X[Drow] += A[(*N)*Dcol - ((Dcol*(Dcol+1))/2) + Drow] * B[Dcol];
    }
  }
   
}

void cor_diag(Real *xloc, Real *yloc, Integer *nl, Integer *cornr, 
	      Real *phi, Real *kappa, Real *res) 
     
     /* This function computes the lower triangle of correlation or distance 
	matrix for a set of  locations.
	(including the diagonal term of the matrix).
	
	cornr defines the correlation function

	0: compute distances only
	otherwise uses corrfctvalue
        xloc, yloc     : xy coordinates of the locations
	nl,            : number of locations
	cornr          : a number indicating the correlation model
	phi            : parameter of the correlation function (scale parameter)
	kappa          : extra parameter for some correlation functions (shape parameter)
	res            : stores the results to be returned, a vector with the 
	lower triangle of the correlation or distance matrix
     */   
     
{ 
  
  Integer register i,j, ind;
  Real register dx, dy;
  Real h;
  
  ind = 0;
  
  for (j=0; j<*nl; j++) {  
    for (i=j; i<*nl; i++) {
      if(i == j){
	if(*cornr > 0) res[ind] = 1.0 ;
	else res[ind] = 0.0 ;
      }
      else{
	dx = (xloc[j] - xloc[i]) ;
	dy = (yloc[j] - yloc[i]) ;
	h  = pythag(dx,dy) ;
	if(*cornr > 0){
	  if(*phi > 0){
	    res[ind] = corrfctvalue((*phi), (*kappa), h, (*cornr));
	  }
	  else res[ind] = 0;
	}
	else
	  res[ind] = h ;
      }
      ind++ ;
    }
  }
}


void kb_sim(Real *means, Real *nscores,
	    Real *lowerA, Real *diagA, 
	    Real *Xmatrix, Integer *Nbig, 
	    Integer *Nsmall, Real *Dval, 
	    Integer *Nsims, Real *Vsqglchi, Real *ss,
	    Real *Blower, Real *Bdiag, Real *bvec, 
	    Integer *Bsize, Real *R0lower) 
{
  /*  
      on exit, for each simulation, nscores = means + chol(varmatrix) %*% nscores
      
      varmatrix was replaced by R0lower in the input      
      Integer Nvarmatrix = *Nbig * (*Nbig +1)/2;
      Real *varmatrix = (Real*) malloc(sizeof(Real)*Nvarmatrix);
  */
  
  lower_R0minusXAXplusBvar(lowerA, diagA, Xmatrix, *Nbig, *Nsmall, Dval, 
			   Blower, Bdiag, bvec, *Bsize, ss, R0lower);  
  
  mvnorm(means, R0lower, nscores, *Nbig, *Nsims, Vsqglchi);
  
  /*  free(varmatrix); */  
}

void mult_kb_sim(Real *means, Real *nscores, Real *lowerA, 
		 Real *diagA, Real *Xmatrix, Integer *Nbig, 
		 Integer *Nsmall, Real *Dval,  Integer *Nsims, 
		 Real *Vsqglchi, Real *ss, Real *Blower, Real *Bdiag, 
		 Real *bvec, Integer *Bsize, Real *R0lower) 
{
  /*    
	OFC: version of kb_sim where means is a matrix 
  */
  
  lower_R0minusXAXplusBvar(lowerA, diagA, Xmatrix, *Nbig, *Nsmall, Dval, 
			   Blower, Bdiag, bvec, *Bsize, ss, R0lower);

  multmvnorm(means, R0lower, nscores, *Nbig, *Nsims, Vsqglchi);  
}

void kb_sim_new(Real *means, Real *nscores,
		Real *lowerA, Real *diagA, 
		Real *Xmatrix, Integer *Npred, 
		Integer *Ndata, Real *Dval, 
		Integer *Nsims, Real *Vsqglchi, Real *ss,
		Real *Blower, Real *Bdiag, Real *bvec, 
		Integer *Bsize, Real *xlocpred, Real *ylocpred,
		Integer *cornr, Real *phi, Real *kappa, 
		Integer *diffmean) 
{
  Integer NR0lower = *Npred * (*Npred + 1)/2 ;
  Real *R0lower = (Real*) malloc(sizeof(Real)*NR0lower) ;
  
  cor_diag(xlocpred, ylocpred, Npred, cornr, phi, kappa, R0lower) ;  
  
  lower_R0minusXAXplusBvar(lowerA, diagA, Xmatrix, *Npred, *Ndata, Dval, 
			   Blower, Bdiag, bvec, *Bsize, ss, R0lower);
  
  if(*diffmean == 0)
    mvnorm(means, R0lower, nscores, *Npred, *Nsims, Vsqglchi);
  else 
    multmvnorm(means, R0lower, nscores, *Npred, *Nsims, Vsqglchi);
  
  free(R0lower);
}


void fastfourier(data,nn,ndim,isign)
     /* 
	taken from the book
	Numerical Recipes in Pascal
	W.H. Press, B.P. Flannery, S.A, Teukolsky, W.T.Vetterling
	Cambridge 1989
     */
     
     double *data;
     int    *nn,ndim,isign;
{ 
  Integer i1,i2,i3,i2rev,i3rev,ibit,idim ;
  Integer  ip1,ip2,ip3,ifp1,ifp2,k1,k2,n ;
  Integer ii1,ii2,ii3 ;
  Integer nprev,nrem,ntot ;
  Integer ip2divip1,ip1div2,ip3divip2,ifp1divip1,Xip1div2,ip3divifp2 ; 
  Real tempi,tempr,wrs,wis,theta,wi,wpi,wpr,wr,wtemp ;
  
  ntot=1; for(idim=0;idim<ndim;idim++){ntot*=nn[idim];} /**/
  nprev=1;
  for (idim=ndim-1;idim>=0;idim--){/**/
    n=nn[idim];
    nrem = ntot / (n*nprev);
    ip1=2 *nprev;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev =1;
    /* */
    ip2divip1 = (ip2-1)/ ip1;
    for (ii2=0;ii2<=ip2divip1;ii2++){
      i2=1+ii2*ip1;
      if (i2<i2rev) {
	ip1div2 = (ip1-2) / 2;
	for (ii1=0; ii1<=ip1div2;ii1++){
	  i1=i2+ii1*2;
	  ip3divip2 = (ip3-i1) / ip2;
	  for (ii3=0;ii3<=ip3divip2;ii3++){
	    i3=i1+ii3*ip2;
	    i3rev = i2rev+i3-i2;
	    tempr = data[i3-1];
	    tempi = data[i3];
	    data[i3-1]=data[i3rev-1];
	    data[i3]=data[i3rev];
	    data[i3rev-1]=tempr;
	    data[i3rev]=tempi;
	  }
	}
      } /* if */
      ibit =ip2 /2;
      while ((ibit >= ip1) &&  (i2rev>ibit)){
	i2rev-= ibit;
	ibit /= 2;
      }
      i2rev+=ibit;
    } /* for */
    /* */
    ifp1 = ip1;
    while (ifp1<ip2) {
      ifp2 = 2*ifp1;
      theta= isign * M_2PI / (ifp2 / ip1);
      wpr  = sin(0.5*theta); wpr*=wpr; wpr*=-2.0;
      wpi  = sin(theta);
      wr   = 1.0;
      wi   = 0.0;
      ifp1divip1 = (ifp1-1) / ip1;
      for (ii3=0;ii3<=ifp1divip1;ii3++) {
	i3 = 1 + ii3*ip1;
	wrs=wr;
	wis=wi;
	Xip1div2=(ip1-2) / 2;
	for (ii1=0;ii1<=Xip1div2;ii1++){
	  i1 = i3+ii1*2;
	  ip3divifp2 = (ip3-i1) / ifp2;
	  for (ii2=0; ii2<=ip3divifp2; ii2++){
	    i2 = i1+ii2*ifp2;
	    k1 = i2;
	    k2 = k1+ifp1;
	    tempr = wrs*data[k2-1]-wis*data[k2];
	    tempi = wrs*data[k2]+wis*data[k2-1];
	    data[k2-1] = data[k1-1]-tempr;
	    data[k2] = data[k1]-tempi;
	    data[k1-1] += tempr;
	    data[k1] += tempi;
	  }
	}
	wtemp = wr;
	wr = wr*wpr-wi*wpi+wr;
	wi = wi*wpr+wtemp*wpi+wi;
      } /* ii3 */
      ifp1 = ifp2;
    } /* while */
    nprev *= n;
  } /* for idim */
}	            

/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper by Wood & Chan 1994)     */
/*********************************************************************/
void woodandchan(cornr, nn, ln, step, phi, kappa, res, m)
     /* 
	cornr : the implemented covariance functions have numbers, see below
	nn    : vector of the number of points in each direction
	ln    : length(nn) == dimension of the randon field
	step  : simulation step (the same in each direction)
	phi   : scale parameter of the correlation finction
	kappa : shape [parameter of the correlation function
	res   : OUT: result vector (tensor)
	m     : IN/OUT: the size of the actually used, enlarged, matrix, 
	see below; IN: if m=(0,...,0) then the matxrix size is 
	automatically determined.
	implemented here only for rotationsinvariant covariance functions
	for arbitrary dimensions;
	(so it works only for even covariance functions in the sense of 
	Wood and Chan,p. 415, although they have suggested a more general 
	algorithm;
	the approximating part has not been implemented either.) 
	If an error appears: m[0]=0 will be returned.
	
     */
     double *step,*res,*phi, *kappa;
     int    *cornr,*nn, *ln, *m;
     
{ 
  int dim, i,j, mtot, maxm,*halfm, Ntot, k, kk, first ;
  int positivedefinite,*cumm,*index,halfmtot,zeroORmiddle ;
  int lastnonzeroORmiddle,modtot,m0eq0 ;
  double *c, *invn2, UU, VV, XX, YY, h, invsqrttwo, epsilon, invsqrtmtot;  

  mtot=0;
  c=0;


  dim= *ln;    /* just renaming */
  switch (dim) {
  case 1 : maxm = 67108864; break;
  case 2 : maxm = 8192; break;
  case 3 : maxm = 256; break;
  case 4 : maxm = 64; break;
  default: error("dim (%d) is too large", dim); return;
  }

  invsqrttwo = 1/sqrt(2.0);
  epsilon = 1e-7;
  m0eq0 = m[0]==0;
  
  halfm = (int*) malloc(sizeof(int) * dim);
  cumm = (int*) malloc(sizeof(int) * dim); 
  /* 
     cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an
     index of the vector (the way the matrix is stored) corresponding
     to the matrix                   
  */
  invn2 = (double*) malloc(sizeof(double) * dim);
  index = (int*) malloc(sizeof(int) * dim);
  
  
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */
  
  for (i=0;i<dim;i++){
    if (m0eq0) {m[i] = 1 << (1+(int) ceil(log(((double) nn[i])-1.0)/log(2.0)));}
    halfm[i] = m[i]/2;
    if (m[i]>maxm) { m[0]=0; return; /* error */} 
    invn2[i] = *step * *step; 
    /* 
       These are the nominators in (3.1).
       But here a rectangle nn[0] * step x ... x nn[dim-1] * step
       is used instead of the [0,1]^d cube.
       "*step" is already squared as finally the Euclidean distance 
       has to calculated. 
       Here, the notation is more general than needed (for extension
       to non-isotropic random fields */
  }
  
  positivedefinite = 0;     
  /* 
     Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
     so only the first to rows of (3.9) (without the taking the
     modulus of h in the first row)
     The following variable `index' corresponds to h(l) in the following
     way: index[l]=h[l]        if 0<=h[l]<=m[l]/2
     index[l]=h[l]-m[l]   if m[l]/2+1<=h[l]<=m[l]-1     
     Then h[l]=(index[l]+m[l]) mod m[l] !!
  */
  for(i=0;i<dim;i++){index[i]=1-halfm[i];} 
  
  /* find Fourier transform of matrix; the latter might be enlarged 
     automatically  */
  while (!positivedefinite){
    cumm[0]=1; for(i=0;i<dim-1;i++){cumm[i+1]=cumm[i] * m[i];} 
    mtot=cumm[dim-1] * m[dim-1]; 
    Rprintf("\n mtot=%d ",mtot);
    /* meaning of following variable c, see eq. (3.8) */
    if ((c =(double*) malloc(2*sizeof(double) * mtot)) ==0){m[0]=0;return;}
    for (i=0;i<mtot;i++){ /* fill in c(h) column-wise*/
      j=0; 
      for (k=0;k<dim;k++) {j+=cumm[k] * ((index[k]+m[k]) % m[k]);}
      
      /* here specialisation to rotationinvariant cov.-functions: */
      h=0; /* Euclidean distance for the right hand side of eq (3.8) */
  
      for (k=0;k<dim;k++) {h+=invn2[k] * (double) (index[k] * index[k]);}
      c[2*j]=corrfctvalue(*phi, *kappa, sqrt(h), *cornr); c[2*j+1]=0; 
      /* cf. eq. 3.8; first real part, then imaginary part */
      
      k=0;while((k<dim)&&(++index[k]>halfm[k])) {index[k]=1-halfm[k]; k++;}
    }

    fastfourier(c,m,dim,-1);   

    /* check if positive definite. If not enlarge and restart */
    i=0;
    while ((i<mtot)&&(positivedefinite=(c[2*i]>=0 && fabs(c[2*i+1])<epsilon)))
      {i++;}
    if (!positivedefinite) {
      Rprintf(" nonpos %d  %f %f  ",i,c[2*i],c[2*i+1]);
      free(c);
      for (i=0;i<dim;i++) {
        m[i] <<= 1; 
        if (m[i]>maxm) {m[0]=0;return;} /* the better strategy (approximation)
                                           is not implemented yet */
      }           
    }
  }
  
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/
  halfmtot =  mtot /2;
  for(i=0;i<dim;i++){index[i]=0;}
  modtot = cumm[dim-1] * (halfm[dim-1]+1); 
  
  GetRNGstate();
  
  for (i=0; i<modtot; i++) { 
    /* about half of the loop is clearly unnecessary, 
       as pairs of r.v. are created, see inequality * below; 
       therefore modtot is used instead of mtot;
       as the procedure is complicated, it is easier to let the loop run
       over a bit more values of i than using exactly mtot/2 values */ 
    zeroORmiddle = 1;
    lastnonzeroORmiddle=0; 
    /* theoretically, it is undefined, if
       index[k] = (0 or halfm[k])  for all k;
       but this would lead an error in the inequality *
       below; in the case all the indices are 0 or
       halfm[k] then inequality * below is obviously
       satisfied with lastnonzeroORmiddle=0
    */            
    kk=0; 
    /* kk takes the role of k[l] in Prop. 3; but here kk is already 
       transformed from matrix indices to a vector index */  
    for (k=0; k<dim; k++){
      if (index[k]==0 || index[k]==halfm[k]) {kk += cumm[k] * index[k];}
      else {
	zeroORmiddle=0;
	kk += cumm[k] * (m[k]-index[k]);
	lastnonzeroORmiddle=k;
      }
    }
    if (index[lastnonzeroORmiddle]<=halfm[lastnonzeroORmiddle]){ /* * */
      /* By Prop. 3, pairs of r.v. exist which are (neg.) dependent;
	 so, these pairs are created together what is done in
	 the following. The above inequality  guaratees, that
	 those pairs are created only once ---
	 
	 Assume m=4 for all k and dimension=3; then the above inequality 
	 becomes clear considering the following pair of tripels whose 
	 corresponding random variables are simulated simultaneously.
	 [(0, 0, 0); (0, 0, 2)];
	 [(0, 1, 0); (0, 3, 0)]; 
	 [(0, 3, 0); ...] is ruled out ! 
	 [(1, 1, 1); (3, 3, 3)];
	 [(1, 3, 1); (3, 1, 3)];
	 [(3, 1, 3); ...] is ruled out!
	 [(3, 3, 3); ...] is ruled out!
	 
	 [(0, 0, 2); ...] is NOT ruled out by *, but by the following 
	 inequality (i<halfmtot)
      */
      
      if (zeroORmiddle){ /* step 2 of Prop. 3 */
        if (i<halfmtot) { 
	  /* create two Gaussian random variables XX & YY  with variance 1/2 */ 
          UU = sqrt(-log(unif_rand()));
          VV = unif_rand() * M_2PI;
          XX = UU * sin(VV);
          YY = UU * cos(VV);
          
          kk = i + halfmtot;
          c[2*i] = sqrt(c[2*i]) * XX; c[2*i+1]=0.0;
          c[2*kk] = sqrt(c[2*kk]) * YY; c[2*kk+1]=0.0;
	}
      } else { /* step 3 of Prop. 3 */
        UU = sqrt(-log(unif_rand()));
        VV = unif_rand() * M_2PI;
        XX = UU * sin(VV);
        YY = UU * cos(VV);
        c[2*i+1]=c[2*i]=sqrt(c[2*i]); c[2*i]*=XX; c[2*i+1]*=YY;
        c[2*kk+1]=c[2*kk]=sqrt(c[2*kk]); c[2*kk]*=XX; c[2*kk+1]*= -YY;
      }
    }
    k=0; while((k<dim) && (++index[k]>=m[k])) {index[k++]=0;}
  }
  
  PutRNGstate();

  fastfourier(c,m,dim,-1);
  
  /* now we correct the result of the fastfourier transformation
     by the factor 1/sqrt(mtot) and read the relevant matrix out of 
     the large vector c */
  Ntot=1; for(i=0;i<dim;i++){Ntot *= nn[i];} 
  invsqrtmtot = 1/sqrt(mtot);
  first = 1;
  for(i=0;i<dim;i++){index[i]=0;}
  for (i=0;i<Ntot;i++){
    kk=0;for (k=0; k<dim; k++){kk+=cumm[k] * index[k];}
    res[i]=c[2*kk]*invsqrtmtot;
    if ((fabs(c[2*kk+1])>epsilon) && (first)){
      Rprintf("IMAGINARY PART <> 0"); first=0;}
    k=0; while((k<dim) && (++index[k]>=nn[k])) {index[k++]=0;}
  }
  free(halfm); free(cumm); free(invn2); free(index); free(c);
}


 

