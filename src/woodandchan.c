#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

#define TWOPI 6.28318530717959


void fastfourier(data,nn,ndim,isign)
   /* taken from the book
      Numerical Recipes in Pascal
      W.H. Press, B.P. Flannery, S.A, Teukolsky, W.T.Vetterling
      Cambridge 1989
   */

   double *data;
   int    *nn,ndim,isign;
{ 
   int i1,i2,i3,i2rev,i3rev,ibit,idim,
       ip1,ip2,ip3,ifp1,ifp2,k1,k2,n,
       ii1,ii2,ii3,
       nprev,nrem,ntot,
       ip2divip1,ip1div2,ip3divip2,ifp1divip1,Xip1div2,ip3divifp2;
   double tempi,tempr,wrs,wis,theta,wi,wpi,wpr,wr,wtemp;

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
	theta= isign * TWOPI / (ifp2 / ip1);
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

/**************************************************************/
/*                  AUXILIARY FUNCTIONS                       */
/**************************************************************/
double besselI0(double x)/* Press&Flannery&Teukolsky&Vetterling, pp 197ff */
{ double ax,y;
  if ((ax = fabs(x)) < 3.75) {    y = x/3.75;    y *=y;
    return (1.0 + y*(3.5156229 + y*(3.0899424 +  y*(1.2067492 + 
		 y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2))))));

  } else{ y = 3.75/ax;
    return ((exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+ 
	   y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+ 
	   y*(-0.1647633e-1+y*0.392377e-2)))))))));
  }
}
double besselI1(double x)/* Press&Flannery&Teukolsky&Vetterling, pp 197ff */
{ double ax,y,ans;
  if ((ax = fabs(x)) < 3.75) { y = x / 3.75; y *= y;
    ans = ax*(0.5 + y*(0.87890594 +  y*(0.51498869 + y*(0.15084934 + 
	   y*(0.2658733e-1 + y*(0.301532e-2 +  y*0.32411e-3))))));
  } else{ y = 3.75/ax;
    ans = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans = 0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+ 
	  y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax) / sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}
double besselK0(double x) /* Press&Flannery&Teukolsky&Vetterling, pp 197ff */
{ double y;
  if (x <= 2.0) { y = x*x/4.0;
    return((-log(x/2.0)*besselI0(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+
          y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5)))))));
  } else{ y = 2.0/x;
    return((exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+ 
	  y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))))); 
  }
}
double besselK1(double x)/* Press&Flannery&Teukolsky&Vetterling, pp 197ff */
{ double y;
  if (x <= 2.0) { y = x * x / 4.0;
     return((log(x/2.0)*besselI1(x))+(1.0/x)*(1.0+y*(0.15443144+
             y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+
             y*(-0.110404e-2+y*(-0.4686e-4))))))));
  } else{ y = 2.0/x;
     return ((exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+ 
	y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3))))))));
  }
}
double besselK(long n, double x) /* Press et al, pp 197ff */
/* returns the modified Bessel function K_n(x) for positive x, n >= 2 */
{ long j;  double bk, bkm, bkp, tox;
  if (n < 0) error("Index nu less than 0 in besselK"); 
  else if (n==0) {return besselK0(x);} 
  else if (n==1) {return besselK1(x);}
  tox = 2.0/x;  bkm = besselK0(x);   bk  = besselK1(x);
  for(j=1; j<n; j++){  bkp = bkm + j * tox * bk;   bkm = bk;  bk  = bkp;}
  return bk;
}

void besselKproc(n,x,K)
   int *n; double *x,*K;
{
   long nn;
   nn = (long) *n;
   K[0] = besselK(nn,*x);
}
     


/**************************************************************/
/*                  COVARIANCE FUNCTIONS                      */
/**************************************************************/

double *PARAM; /* Parameter for the covariance functions
                  the first is the scale, 
                  the second the  shape parameter,
                  then further parameters */
                  

/* correlation functions : */
double nugget(double x){return (x==0);}

double spherical(double x){
  double y;
  if (x>PARAM[0]) {return 0.0;}
  else {y=fabs(x/PARAM[0]); 
        return(1.0+y*0.5*(y*y-3.0));}
}

double stable(double x){
  if (x==0) {return(1.0);}
  else {return(exp(-1.0*exp(log(fabs(x/PARAM[0]))*PARAM[1])));}
}
  

double whittlematern(double x)
{ long i,nu; double ans,y;
  if (x==0) {return 1;}
  y = fabs(x/PARAM[0]);
  nu = (long) PARAM[1];
  ans = y*besselK(nu,y);
  for (i=1;i<nu;i++){ans*=y/2/i;}
  return ans;
}

void Xwhittlematern(x,nx,phi,nu)
  double *x,*phi,*nu;
  long *nx;
{ long i,j; double ans,y;
  for (j=0;j<*nx;j++){
    if (x[j]==0) {x[j]=1;} else {
      y = fabs(x[j]/ *phi);
      ans = y*besselK(*nu,y);
      for (i=1;i< *nu;i++){ans*=y/2/i;}
      x[j]=ans;
    } 
  }
}
  
  

double gneiting(double x){ 
  double y,y2,oneMy8;
  if (x>PARAM[0]) {return(0.0);}
  else {y = fabs(x/PARAM[0]);
        y2= y*y;
        oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
        return ((1.0+8.0*y+25.0*y2+32.0*y2*y)*oneMy8);}
}


double gneitingdiff(double x){ 
  /* three Parameters!  
     overall scale 
     whittle-matern shape
     gneiting scale
  */
  double y,y2,oneMy8;
  y = PARAM[2]*PARAM[0];
  if (x>y) {return(0.0);}
  else {y = fabs(x/y);
        y2= y*y;
        oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
        return ((1.0+8.0*y+25.0*y2+32.0*y2*y)*oneMy8*whittlematern(x)); }
}


double cauchy(double x){
  double y;
  y = x/PARAM[0];
  return pow(1+y*y,-PARAM[1]);
}

double holeeffect(double x){
  double y;
  if (x==0) {return 1;}
  y = x/PARAM[0];
  return sin(y)/y;
}

/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/
void woodandchan(covnr,nn,ln,step,param,lparam,res, m)
  /* covnr : the implemented covariance functions have numbers, see below
     nn    : vector of the number of points in each direction
     ln    : length(nn) == dimension of the randon field
     step  : simulation step (the same in each direction)
     param : most of the covariance function have parameters, 
             starting with scale parameter, then shape parameter (check 
             definitions above for more complicated models like Gneiting's 
             family!)
     lparam: number of parameters
     res   : OUT: result vector (tensor)
     m     : IN/OUT: the size of the actually used, enlarged, matrix, 
             see below; IN: if m=(0,...,0) then the matrix size is 
             automatically determined.
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;
     the approximating part has not been implemented either.) 
     If an error appears: m[0]=0 will be returned.
  */
  double *step,*res,*param;
  int    *covnr,*nn, *ln, *m, *lparam;
{ 
  int     dim, i,j, mtot, maxm,*halfm, Ntot, k, kk, first,
          positivedefinite,*cumm,*index,halfmtot,zeroORmiddle,
          lastnonzeroORmiddle,modtot,m0eq0;
  double *c,*invn2,UU,VV,XX,YY,h,invsqrttwo, epsilon,invsqrtmtot;
  double (*cov)(double);

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
 
  PARAM = (double*) malloc(sizeof(double) * *lparam);
  for(i=0;i< *lparam; i++) {PARAM[i] = param[i];}
  switch (*covnr) {
  case 1 : cov=spherical; break;
  case 2 : cov=gneiting; break;
  case 3 : cov=stable; break;
  case 4 : cov=whittlematern; break;
  case 5 : cov=gneitingdiff; break;
  case 6 : cov=cauchy; break;
  case 7 : cov=holeeffect; break;
  default: cov=nugget;
  }
 
  halfm = (int*) malloc(sizeof(int) * dim);
  cumm = (int*) malloc(sizeof(int) * dim); /* cumm[i+1]=\prod_{j=0}^i m[j] 
          cumm is used for fast transforming the matrix indices into an
          index of the vector (the way the matrix is stored) corresponding
          to the matrix                   */
  invn2 = (double*) malloc(sizeof(double) * dim);
  index = (int*) malloc(sizeof(int) * dim);


  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */
  
  for (i=0;i<dim;i++){
    if (m0eq0) {m[i] = 1 << (1+(int) ceil(log(((double) nn[i])-1.0)/log(2.0)));}
    halfm[i] = m[i]/2;
    if (m[i]>maxm) { m[0]=0; return; /* error */} 
    invn2[i] = *step * *step; /* These are the nominators in (3.1).
         But here a rectangle nn[0] * step x ... x nn[dim-1] * step
         is used instead of the [0,1]^d cube.
         "*step" is already squared as finally the Euclidean distance 
         has to calculated. 
         Here, the notation is more general than needed (for extension
         to non-isotropic random fields */
  }

  positivedefinite = 0;     
  /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
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
      c[2*j]=cov(sqrt(h)); c[2*j+1]=0; /* cf. eq. 3.8; first real part, 
                                          then imaginary part */

      k=0;while((k<dim)&&(++index[k]>halfm[k])) {index[k]=1-halfm[k]; k++;}
    }

    fastfourier(c,m,dim,-1);   
    /* check if possitive definite. If not enlarge and restart */
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

  for (i=0; i<modtot; i++) { /* about half of the loop is clearly unnecessary, 
          as pairs of r.v. are created, see inequality * below; 
          therefore modtot is used instead of mtot;
          as the procedure is complicated, it is easier to let the loop run
          over a bit more values of i than using exactly mtot/2 values */ 
    zeroORmiddle = 1;
    lastnonzeroORmiddle=0; /* theoretically, it is undefined, if
                              index[k] = (0 or halfm[k])  for all k;
                              but this would lead an error in the inequality *
                              below; in the case all the indices are 0 or
                              halfm[k] then inequality * below is obviously
                              satisfied with lastnonzeroORmiddle=0
                           */            
    kk=0; /* kk takes the role of k[l] in Prop. 3; but here kk is already 
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
          VV = unif_rand() * TWOPI;
          XX = UU * sin(VV);
          YY = UU * cos(VV);
          
          kk = i + halfmtot;
          c[2*i] = sqrt(c[2*i]) * XX; c[2*i+1]=0.0;
          c[2*kk] = sqrt(c[2*kk]) * YY; c[2*kk+1]=0.0;
	}
      } else { /* step 3 of Prop. 3 */
        UU = sqrt(-log(unif_rand()));
        VV = unif_rand() * TWOPI;
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
  free(PARAM);
}

 

