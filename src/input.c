#include "utilities.h"

void malloc_gene_data(GENE_DATA *pdata)
{
  int i;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  
  assert(pdata->d=(float **)malloc(sizeof(float*)*nrow));
  assert(pdata->L=(int*)malloc(sizeof(int)*ncol));

  /*initialization*/
  memset(pdata->L,0,sizeof(int)*ncol);
  for(i=0;i<ncol;i++) 
    pdata->L[i]=0;

  for (i=0; i<nrow; i++) {
    assert(pdata->d[i]=(float *) malloc(sizeof(float)*ncol));
  }
}


/*************************************************************************/
/*                           free_gene_data                              */
/*************************************************************************/
/*free the space allocated for pdata*/
void free_gene_data(GENE_DATA *pdata)
{
  int i;
  for (i=0; i<pdata->nrow; i++) {
    free(pdata->d[i]);
  }
  free(pdata->L);
  free(pdata->d);
}


void create_gene_data(double*d,  int*pnrow, int*pncol, int*L, GENE_DATA* pdata)
{
  int i,j, nL=0;
  pdata->nrow=*pnrow;
  pdata->ncol=*pncol;
  malloc_gene_data(pdata);
  for (j=0; j<pdata->ncol; j++) 
    pdata->L[j]=L[j];

  for (i=0;i<*pncol;i++) {
    if (L[i]>nL)
      nL=L[i];
  }
  nL++;
  pdata->nL=nL;
  
  for (i=0; i<pdata->nrow; i++) {
    for (j=0; j<pdata->ncol; j++) {
      pdata->d[i][j]=d[j*pdata->nrow+i];
    }
  }
}


/********************************************************************************/
/*             print_gene_data                                                  */
/********************************************************************************/
/*print the gene_data to stderr, useful in debug*/
void print_gene_data(GENE_DATA *pdata)
{
  int i, j;
  for (i=0; i<pdata->nrow; i++){
    for (j=0; j<pdata->ncol; j++) 
      Rprintf(" %5.3f", pdata->d[i][j]);
    Rprintf("\n");
  }
}

void create_deds_res(int *pnrow, int *pnsig, int *pnT, DEDS_RES *pdr)
{
  int i;
  pdr->nrow=*pnrow;
  pdr->nsig=*pnsig;
  pdr->nT=*pnT;

  assert(pdr->R=(int *)malloc(sizeof(int)*(*pnrow)));
  assert(pdr->E=(float *)malloc(sizeof(float)*(*pnT)));
  assert(pdr->D=(float *)malloc(sizeof(float)*(*pnrow)));
  assert(pdr->wval=(float *)malloc(sizeof(float)*(*pnT)));
  assert(pdr->FDR=(double *)malloc(sizeof(double)*(*pnrow)));
  assert(pdr->T=(float **)malloc(sizeof(float *)*(*pnrow)));
  for(i=0;i<(*pnrow);i++)
    assert(pdr->T[i]=(float *)malloc(sizeof(float)*(*pnT)));
  
}
 
void free_deds_res(DEDS_RES *pdr)
{
  int i;
  free(pdr->R);
  free(pdr->E);
  free(pdr->D);
  free(pdr->FDR);
  for(i=0;i<(pdr->nrow);i++)
    free(pdr->T[i]);
  free(pdr->T);
  
}
 
void create_tmod_data(int *pnrow, TMOD_DATA *ptmod)
{
  ptmod->nrow=*pnrow;
  
  assert(ptmod->mean=(float *)malloc(sizeof(float)*(*pnrow)));
  assert(ptmod->sigma2=(float *)malloc(sizeof(float)*(*pnrow)));
  assert(ptmod->df_resid=(int *)malloc(sizeof(int)*(*pnrow)));
  assert(ptmod->stdev_unscale=(float *)malloc(sizeof(float)*(*pnrow)));
    
}
 
void free_tmod_data(TMOD_DATA *ptmod)
{
  free(ptmod->mean);
  free(ptmod->sigma2);
  free(ptmod->df_resid);
  free(ptmod->stdev_unscale);
}
 
/********************************************************************************/
/*                    sort_gene_data                                            */
/********************************************************************************/
/* Description:
   sort the rows of gene_data such that row R[i] of is the first row, i=0,...,m-1,
   wher R[0],...,R[m-1] is a permutation of (0,...,m-1)
*/
										 
void sort_gene_data(GENE_DATA* pdata,int*R)
{
  int i,nrow=pdata->nrow;
  float** old_d;  /*th old addresses of the gene data*/
  assert(old_d=(float**)malloc(sizeof(float*)*nrow));
  /*store the original pointers from pdata*/ 
  for(i=0;i<nrow;i++)
    {
      old_d[i]=pdata->d[i];
    }
  /*rearrange the data so that it's ordered according to R*/
  for(i=0;i<nrow;i++)
    {
      pdata->d[i]=old_d[R[i]];
    }
  free(old_d);
}
/********************************************************************************/
/*                    sort_vector                                               */
/********************************************************************************/
/* Desciption
      sort the vector V according to the order R with n elemnets
      where R[0],...,R[n-1] is a permutation of 0,...n-1
*/
void sort_vector(float* V,int*R,int n)
{
  float* old_V;
  int i;
  assert(old_V=(float*)malloc(sizeof(float)*n));
  for(i=0;i<n;i++)
    old_V[i]=V[i];
  for(i=0;i<n;i++)
    V[i]=old_V[R[i]];
  free(old_V);
}
