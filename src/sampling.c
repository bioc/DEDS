#include "utilities.h"


static int l_n=0; /* the length of L */
static int l_k=0; /* the number of groups */
static int *l_nk=NULL; /* the number of objects in each group */
static int *l_L=NULL; /* first label */
static int l_b=0; /* the number of permutations done so far */
static int l_B=0; /*  the number of all permuations */
static int *l_permun=NULL;
static int *l_ordern=NULL;
/*extern long int g_random_seed;*/
long int g_random_seed=3455660;

void creat_sampling(int n, int *L, int B)
{
  int i,k;
  l_n=n;
  l_B=B;
  l_b=0;
  
  /* initialize l_L which stores the original label L */
  assert(l_L=(int *)malloc(n*sizeof(int)));
  memcpy(l_L, L, sizeof(int)*n);

  /* calculate number of classes */
  k=0;
  for (i=0;i<n;i++) {
    if (L[i]>k)
      k=L[i];
  }
  k++;
  l_k=k;
  
  /* calculate number of cases in each class */
  assert(l_nk=(int *)malloc(k*sizeof(int)));
  memset(l_nk, 0, sizeof(int)*k);
  for (i=0;i<n;i++){
    l_nk[L[i]]++;
  }
      
  /* setting the maximum B*/
  /* this is checked in R already, we don't need it here */
  /*maxB=1;
  rest=n;
  if (l_k>1) {
    for(i=0;i<l_k;i++){
      maxB*=bincoeff(rest,l_nk[i]);
      rest-=l_nk[i];
    }
  }
  else 
    maxB=(int)dpowern(2.0, l_n);
    
  
  if(B>=maxB) {
    l_B=maxB;
    Rprintf("\nWe are doing %d complete permutations\n", l_B);
  }
  else {
    l_B=B;
    Rprintf("\nWe are doing %d random permutations\n", l_B);
    }*/

  /* initialize l_permun and l_ordern */
  assert(l_permun=(int *)malloc(n*sizeof(int)));
  assert(l_ordern=(int *)malloc(n*sizeof(int)));
  for (i=0;i<n;i++)
    l_ordern[i]=i;

  set_seed(g_random_seed);
}

void delete_sampling()
{
  free(l_L);
  l_L=NULL;
  free(l_nk);
  l_nk=NULL;
  free(l_permun);
  l_permun=NULL;
  free(l_ordern);
  l_ordern=NULL;
}

int next_sample(int *L)
{
  int n=l_n;
  if (l_b>=l_B) return 0;
  memcpy(l_permun, l_ordern, sizeof(int)*n);
  sample(l_permun, n, n);
  sample2label(n, l_k, l_nk, l_permun, L);
  l_b++;
 
  return(1);
}

int next_sample_1(int *L)
{
  int n=l_n, i;
  if (l_b>=l_B) return 0;
  memcpy(l_permun, l_ordern, sizeof(int)*n);
  for (i=0;i<n;i++) {
    if (get_rand()<0.5) L[i]=1;
    else L[i]=-1;
  }
  l_b++;

  return(1);
}


void sample2label(int n, int k, int *nk,int *permun, int *L)
{
  int l,s,j;
  s=0;/*s is for starting*/
  for(l=0;l<k;l++){
    for(j=0;j<nk[l];j++){
      L[permun[s]]=l;
      s++;
    }
  }
 
} 
      
    

  

   

