#include "utilities.h"

static float* gp_arr;   
/***********************************************************************/
void order_data(float* V,int*R,int n,FUNC_CMP func_cmp)
{
  int i;
  for(i=0;i<n;i++)
    R[i]=i;
  gp_arr=V;
  qsort(R,n,sizeof(R[0]),func_cmp);
}
  

/********************************************************************************/
/*               comparing functions                                            */
/********************************************************************************/
/*  1 cmp_abs: comparing the absolute values, 
    2 cmp_low: comparing the lower tail
    3 cmp_high: comparing the higher tail
  *Note the gp_arr is the global pointer which has to be used in the program qsort
    After using the quick with those functions, it will be odered such that
    1 in cmp_abs: the bigger values in absolute values will in lower index
                  similar to two sided-test
    2 in cmp_low: the smaller values will in lower index
                  similar to lower tail test
    3 in cmp_high:the bigger values will in lower index 
                  similar to hight tail test.
*/
/*always put the absolute value at the end of the array*/
int cmp_abs(const void *v1, const void *v2) {
  float f1=fabs(*(gp_arr+*(int *)v1));
  float f2=fabs(*(gp_arr+*(int *)v2));
  if(f1==NA_FLOAT)
    return 1;
  if(f2==NA_FLOAT)
    return -1;
  if (f1<f2)
    return 1;
  if (f1>f2)
    return -1;
  else
    return 0; 
}
int cmp_low(const void *v1, const void *v2) {
  if((*(gp_arr+*(int *)v1))==NA_FLOAT)
    return 1;
  if((*(gp_arr+*(int *)v2))==NA_FLOAT)
    return -1;
  if ((*(gp_arr+*(int *)v1))<(*(gp_arr+*(int *)v2)))
    return -1;
  if ((*(gp_arr+*(int *)v1))>(*(gp_arr+*(int *)v2)))
    return 1;
  else
    return 0; 
} 
int cmp_high(const void *v1, const void *v2) {
  if((*(gp_arr+*(int *)v1))==NA_FLOAT)
    return -1;
  if((*(gp_arr+*(int *)v2))==NA_FLOAT)
    return 1;
  if ((*(gp_arr+*(int *)v1))<(*(gp_arr+*(int *)v2)))
    return 1;
  if ((*(gp_arr+*(int *)v1))>(*(gp_arr+*(int *)v2)))
    return -1;
  else
    return 0; 
}


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
  dX=(int *)malloc(n*sizeof(int));
  
  for(i=0;i<n;i++) {
    if(R_FINITE(X[i])) {
      dX[total]=i;
      total+=1;
   }
  }
  k=total/2;

  sX=(float *)malloc(total*sizeof(float));
  
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
   dX=(int *)malloc(nX*sizeof(int));
   index=(double *)malloc(nq*sizeof(double));
   lo=(double *)malloc(nq*sizeof(double));
   hi=(double *)malloc(nq*sizeof(double));
   qX=(float *)malloc(nq*sizeof(float));
  
   /* total is number of finite data */
   for(i=0;i<nX;i++) {
    if(R_FINITE(X[i])) {
      dX[total]=i;
      total++;
    }
   }
  
   sX=(float *)malloc(total*sizeof(float));
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

