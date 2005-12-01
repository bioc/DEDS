#include <R_ext/Random.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* codes modified from the multtest package */
/*get the n samples from the n-dim vector V. the results are stored 
in the first m member of vector V*/
void sample(int *V, int n, int m)
{
  GetRNGstate();
  int i,j,temp;
  float f;
  for(i=0;i<m;i++){
    /* no need to worry yet    
       if(i==(n-1)) continue; */ /*no need to swap with the last elements*/
    j=n;
    while (j==n){/*skip the border, we only want random
		    numbers from i,i+1,i+2,...,n-1*/
      f=unif_rand()*(n-i);
      j=i+floor(f);
    }
    /*swap the nubmer V[i] and V[j] whther if i==j*/
    temp=V[j];
    V[j]=V[i];
    V[i]=temp;
  }   
  PutRNGstate();
}
