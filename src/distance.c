#include "utilities.h"
float mad(float *X, int n)
{
  int i;
  float med, res, *ad;
  ad=(float *)malloc(sizeof(float)*n);

  med=median(X, n);
    
  for (i=0;i<n;i++) {
      if (!R_FINITE(X[i])) ad[i]=NA_REAL;
      else  ad[i]=fabs(X[i]-med);
  }
  res=1.4826*median(ad, n);
  return(res);
}



    
void compute_euclid(float **X, int nrow, int ncol, float *E, float *wval, float *dist)
{
  float dev, *D;
  int i, j, *count;

  count=(int *)malloc(nrow*sizeof(int));
  memset(count, 0, nrow*sizeof(int));
  D=(float *)malloc(nrow*sizeof(float));
  memset(D, 0, nrow*sizeof(float));
    
  
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(X[i][j]))
	continue;
      dev=X[i][j]-E[j];
      D[i]+=wval[j]*dev*dev;
      count[i]++;
    }
  }
  
  for (i=0;i<nrow;i++) {
    if (count[i]==0) D[i]= NA_REAL;
    else if (count[i]!=ncol) {
      D[i]/=((float)count[i]/ncol);
      D[i]=sqrt(D[i]);
    }
    else D[i]=sqrt(D[i]);
    dist[i]=D[i];
    }
}


static float *gp_arr;
void order_index(float *V,int *R,int n)
{
  int i;
  float *oV;

  oV=(float *)malloc(sizeof(float)*n);

  for(i=0;i<n;i++){
    R[i]=i;
    oV[i]=V[i];
  }
      
  gp_arr=V;
  qsort(R,n,sizeof(R[0]),indexCompare);
 
}
  
void order_dist(float *V,int n)
{
  gp_arr=V;
  qsort(V,n,sizeof(V[0]),distCompare);
  
}
  
int indexCompare(const void *v1, const void *v2) 
{
  if(!R_FINITE(*(gp_arr+*(int *)v1)))
    return 1;
  if(!R_FINITE(*(gp_arr+*(int *)v2)))
    return -1;
  if ((*(gp_arr+*(int *)v1))<(*(gp_arr+*(int *)v2)))
    return -1;
  if ((*(gp_arr+*(int *)v1))>(*(gp_arr+*(int *)v2)))
    return 1;
  else
    return 0; 
} 

int distCompare(const void *v1, const void *v2) 
{
  const float *i=v1;
  const float *j=v2;
  if(!R_FINITE(*i))
    return 1;
  if(!R_FINITE(*j))
    return -1;
  if ((*i)<(*j))
    return -1;
  if ((*i)>(*j))
    return 1;
  else
    return 0; 
} 
