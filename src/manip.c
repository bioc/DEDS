#include "utilities.h"

/***********************************************************************/
float max_high(float *X, int n)
{
  int i, count=0;
  double tmp=-DBL_MAX;

  for (i=0;i<n;i++){
    if (!R_FINITE(X[i]))
      continue;
    if (X[i]<tmp) 
      continue;
    else 
      tmp=X[i];
    count++;
  }
  
  if (count==0) return(NA_REAL);
  else return(tmp);
}

float max_low(float *X, int n)
{
  int i, count=0;
  double tmp=DBL_MAX;

  for (i=0;i<n;i++){
    if (!R_FINITE(X[i]))
      continue;
    if (X[i]>tmp) 
      continue;
    else 
      tmp=X[i];
    count++;
  }
  if (count==0) return(NA_REAL);
  else return(tmp);
}

float max_abs(float *X, int n)
{
  int i, count=0;
  float tmp=0;

  for (i=0;i<n;i++){
    if (!R_FINITE(X[i]))
      continue;
    if (fabs(X[i])<tmp) 
      continue;
    else 
      tmp=fabs(X[i]);
    count++;
  }
  if (count==0) return(NA_REAL);
  else return(tmp);
}
/********************************************************************************/
/*                      bincoeff                                                */
/********************************************************************************/
/*Descriptions: return the binomial coefficient of n choosing k*/
int bincoeff(int n, int k)
{
  float f=n;
  int i;
  for(i=1;i<k;i++)
    f*=(n-i)/(i+1.0);
  return (int)(f+0.5);
}

float median(float *X, int n)
{
  float *sX, m;
  int *dX;
  unsigned long total=0, i, k;
  assert(dX=(int *)malloc(n*sizeof(int)));
  
  for(i=0;i<n;i++) {
    if(R_FINITE(X[i])) {
      dX[total]=i;
      total+=1;
   }
  }
  k=total/2;

  assert(sX=(float *)malloc(total*sizeof(float)));
  
  for(i=0;i<total;i++)   sX[i]=X[dX[i]];
 
  m=sel(k, total, sX);
  free(sX);
  free(dX);
  return(m);
}

void quantile(float *X, int nX, float *q, int nq, float *ret)
{
   float *sX, *qX;
   double *index, *lo, *hi;
   int *dX;
   unsigned long total=0, i;
   assert(dX=(int *)malloc(nX*sizeof(int)));
   assert(index=(double *)malloc(nq*sizeof(double)));
   assert(lo=(double *)malloc(nq*sizeof(double)));
   assert(hi=(double *)malloc(nq*sizeof(double)));
   assert(qX=(float *)malloc(nq*sizeof(float)));
  
   /* total is number of finite data */
   for(i=0;i<nX;i++) {
    if(R_FINITE(X[i])) {
      dX[total]=i;
      total++;
    }
   }
  
   assert(sX=(float *)malloc(total*sizeof(float)));
   for(i=0;i<total;i++) sX[i]=X[dX[i]];
   qsort(sX, total, sizeof(sX[0]), distCompare);
   
   for(i=0;i<nq;i++) {
     index[i]=(double)(total-1.0)*q[i];
     lo[i]=floor(index[i]);
     hi[i]=ceil(index[i]);
     qX[i]=sX[(int)lo[i]];
   }
   for(i=0;i<nq;i++) {
     if(index[i]==lo[i]) ret[i]=qX[i];
     else ret[i]=qX[i]+(sX[(int)hi[i]]-sX[(int)lo[i]])*(index[i]-lo[i]);
   }

   free(sX);
   free(qX);
   free(index);
   free(lo);
   free(hi);
   free(dX);
}
   

float dpowern(float n, int k)
{
  float f=1.0;
  int i;
  for (i=0;i<k;i++)   f*=n;
  return(f);
}
void test_select(float *X, unsigned long *n)
{
  float x;
  x=median(X, *n);
  Rprintf("%5.3f\n", x);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float sel(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i, ir, j, l, mid;
  float a, temp;
  
  l=0;
  ir=n-1;
  for(;;) {
    if(ir<=l+1) {
      if(ir==l+1 && arr[ir]<arr[l]) {
	SWAP(arr[l],arr[ir])
      }
      return(arr[k]);
    }
    else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid], arr[l+1])
	if(arr[l]>arr[ir]) {
	  SWAP(arr[l], arr[ir])
	}
        if(arr[l+1]>arr[ir]) {
	  SWAP(arr[l+1], arr[ir])
	}
	if(arr[l]>arr[l+1]) {
	  SWAP(arr[l], arr[l+1])
	}
	i=l+1;
	j=ir;
	a=arr[l+1];
	for(;;) {
	  do i++; while(arr[i]<a);
	  do j--; while(arr[j]>a);
	  if(j<i) break;
	  SWAP(arr[i],arr[j])
	}
	arr[l+1]=arr[j];
	arr[j]=a;
	if(j>=k) ir=j-1;
	if(j<=k) l=i;
    }
  }
}

