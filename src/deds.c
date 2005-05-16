#include "utilities.h"

/***********************************************************************************/
/** Name: get_deds_order                                                          **/
/** Purpose: calculate gene order for differential expression given a set         **/
/**          statistics according to the DEDS algorithm                           **/
/** Input: d -- data vector;                                                      **/
/**        pnrow, pncol -- dimensions of data matrix to be formatted into         **/
/**        L -- class labels in 0s and 1s                                         **/
/**        options -- statistics, distance options                                **/
/**        nL -- number of classes; nT -- number of statistics                    **/
/**        B -- number of permutations                                            **/
/**        E -- extreme statistics to be returned                                 **/
/**        R -- gene order to be returned                                         **/
/**        T -- statistics to be returned                                         **/
/***********************************************************************************/

void get_deds_order(double *d, int *pnrow, int *pncol, int *L, char **options, float *extras, 
		   int *nL, int *nT, int *B, double *E, int *R, double *T)
{
  GENE_DATA data;
  TEST_DATA td;
  DEDS_RES dr;
  double *F;

  assert(F=(double *)malloc(sizeof(double)*(*pnrow)));

  create_gene_data(d,pnrow,pncol,L,&data);
  if(type2test(options, &td, nT, nL, extras)==0) return;
  create_deds_res(pnrow, pncol, nT, &dr);

  func_get_order(&data, &td, &dr, B);
  extract_deds_res(&dr, E, R, F, T);
  
  free_gene_data(&data);
}

/******************************************************************************************/
/** Name: get_deds_FDR                                                                   **/
/** Purpose: calculate gene order and FDR for differential expression given a            **/
/**          set of statistics according to the DEDS algorithm                           **/
/** Inputs: similar to get_deds_order, except                                            **/
/**         quick -- if quick, one set permutation is performed for both                 **/
/**                  generating gene order and gene FDR; the permutation dataset         **/
/**                  is stored, calling for a bigger memory. IF not quick,               **/
/**                  permuations are done separately for the generation of gene          **/
/**                  order and gene FDR, therefore no need to store any                  **/
/**                  intermediate results, but it is slower. As a fixed seed             **/
/**                  random number generator is used, results are same for both          **/
/**                  quick or slow.                                                      **/
/**         pnsig -- number of significant genes to be outputed                          **/
/**         FDR -- basically Q values for significant genes assuming all hypotheses null **/
/******************************************************************************************/

void get_deds_FDR(double *d, int *pnrow, int *pncol, int *L, char **options, float *extras, int *quick, 
		   int *nL, int *nT, int *B, int *pnsig, double *E, int *R, double *FDR, double *T)
{
  GENE_DATA data;
  TEST_DATA td;
  DEDS_RES dr;
          
  create_gene_data(d,pnrow,pncol,L,&data);
  if(type2test(options, &td, nT, nL, extras)==0) return;
  create_deds_res(pnrow, pnsig, nT, &dr);
   
  if (*quick) func_deds_quick(&data, &td, &dr, B);
  else {
    func_get_order(&data, &td, &dr, B);
    func_get_FDR(&data, &td, &dr, B);
    }
  
  extract_deds_res(&dr, E, R, FDR, T);
  free_deds_res(&dr);
    free_gene_data(&data);
}


void func_deds_quick(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B) 
{
  int i, j; /* iteration index */
  int b=0, *L, *bL, is_next; /* sampling variables */
  float wval;
  float **T, ***bT, *tmpT;
  float *bE, *E, *fF;
  float *bD, *bMD;
  int ncol=pdata->ncol;
  int nrow=pdata->nrow;
  int nT=pdr->nT;
  int nsig=pdr->nsig;
  int weighted_dist=ptd->weighted_dist;
  FUNC_COMPUTE_STAT func_compute_stat;
  FUNC_COMPUTE_P func_compute_p=ptd->func_compute_p;
  FUNC_MAX func_max=ptd->func_max;
  FUNC_SAMPLE func_next_sample=ptd->func_next_sample;
  float *extras;
  
  assert(extras=(float *)malloc(nT*sizeof(float)));
  memcpy(extras, ptd->extras, nT*sizeof(float));
  assert(bL=(int *)malloc(ncol*sizeof(int)));
  assert(L=(int *)malloc(ncol*sizeof(int)));
  memcpy(L, pdata->L, ncol*sizeof(int));
  assert(T=(float **)malloc(sizeof(float*)*nrow));
  for (i=0;i<nrow;i++)  
    assert(T[i]=(float *)malloc(sizeof(float)*nT));
  assert(tmpT=(float *)malloc(sizeof(float)*nrow));
  assert(bE=(float *)malloc(sizeof(float)*nT));
  assert(E=(float *)malloc(sizeof(float)*nT));
  assert(fF=(float *)malloc(sizeof(float)*nrow));
  assert(bD=(float *)malloc(sizeof(float)*nrow));
  
  assert(bMD=(float *)malloc(sizeof(float)*(nrow*(*B))));
  assert(bT=(float ***)malloc(sizeof(float **)*(*B)));
  for (i=0;i<(*B);i++) {
    assert(bT[i]=(float **)malloc(sizeof(float *)*nrow));
    for (j=0;j<nrow;j++)
      assert(bT[i][j]=(float *)malloc(sizeof(float)*nT));
  }

  /********************************************************/
  /***   compute E by permutation                       ***/
  /********************************************************/

  /* compute the original stat */
  Rprintf("\nE of the orginial data is: "); 
  for (i=0;i<nT;i++) {
    func_compute_stat=ptd->stat_array[i];
    (*func_compute_stat)(pdata, L, tmpT, &extras[i]);
    for (j=0;j<nrow;j++){
      pdr->T[j][i]=tmpT[j];
      if(func_max==max_abs) tmpT[j]=fabs(tmpT[j]);
      T[j][i]=tmpT[j];
    }
    E[i] = (*func_max)(tmpT, nrow);
    Rprintf("%5.3f  ", E[i]);
    if (weighted_dist) {
      wval=mad(tmpT, nrow);
      wval=1/(wval*wval);
      pdr->wval[i]=wval;
    }
    else pdr->wval[i]=1;
    }
  Rprintf("\n");
   
  /* start permutation */
  creat_sampling(ncol, L, (*B));
  is_next=(*func_next_sample)(bL);
  b=0;
  while(is_next){
    for(i=0;i<nT;i++) {
      func_compute_stat = ptd->stat_array[i];
      (*func_compute_stat)(pdata, bL, tmpT, &extras[i]);
      bE[i]=(*func_max)(tmpT, nrow);
      for(j=0;j<nrow;j++) {
	  bT[b][j][i]=tmpT[j];
	  if(func_max==max_abs) bT[b][j][i]=fabs(bT[b][j][i]);
      }
      
      if((func_max==max_high)&&(bE[i]>E[i]))
	E[i]=bE[i];
      else if ((func_max==max_low)&&(bE[i]<E[i]))
	E[i]=bE[i];
      else if ((func_max==max_abs)&&(fabs(bE[i])>fabs(E[i])))
	E[i]=bE[i];
    }
    b++;
    print_b(b,*B,"b=");
    is_next=(*func_next_sample)(bL);
    }

  Rprintf("\nAfter permutation , E is set at: "); 
  for (i=0;i<nT;i++) {
    Rprintf("%5.3f  ", E[i]);
    pdr->E[i]=E[i];
  }
  Rprintf("\n");
    
  compute_euclid(T, nrow, nT, E, pdr->wval, pdr->D); 
  
  order_index(pdr->D, pdr->R, nrow);      
  Rprintf("\nSummarizing DEDS results for %d permutations and %d genes, please wait... \n", (*B), nsig);
	 
  
  for (b=0;b<(*B);b++) {
    compute_euclid(bT[b], nrow, nT, pdr->E, pdr->wval, bD);
    /*qsort(bD, nrow, sizeof(bD[0]), distCompare);*/
    for(i=0;i<nrow;i++) 
      bMD[b*nrow+i]=bD[i];
  }
  
  (*func_compute_p)(bMD, pdr->D, pdr->R, &nrow, B, &nsig, fF);

  for(i=0;i<nrow;i++) pdr->FDR[i]=(double)fF[i];
			    
  
  free(bL);
  free(L);
  free(bE);
  free(E);
  free(fF);
  free(bD);
  free(bMD);
  free(tmpT);
  for (i=0;i<(*B);i++){
    for (j=0;j<nrow;j++)
      free(bT[i][j]);
    free(bT[i]);
  }
  free(bT);
  delete_sampling();
}


void func_get_order(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B)
{
  int i, j; /* iteration index */
  int b=0, *L, *bL, is_next; /* sampling variables */
  float **T, *tmpT;
  float wval;
  float *bE, *E;
  float *D;
  int ncol=pdata->ncol;
  int nrow=pdata->nrow;
  int nT=ptd->n_stat;
  int weighted_dist=ptd->weighted_dist;
  FUNC_COMPUTE_STAT func_compute_stat;
  FUNC_MAX func_max=ptd->func_max;
  FUNC_SAMPLE func_next_sample=ptd->func_next_sample;
  float *extras;
  
  assert(extras=(float *)malloc(nT*sizeof(float)));
  memcpy(extras, ptd->extras, nT*sizeof(float));
  assert(bL=(int *)malloc(ncol*sizeof(int)));
  assert(L=(int *)malloc(ncol*sizeof(int)));
  memcpy(L, pdata->L, ncol*sizeof(int));
  assert(tmpT=(float *)malloc(sizeof(float)*nrow));
  assert(T=(float **)malloc(sizeof(float *)*nrow));
  for (i=0;i<nrow;i++)
    assert(T[i]=(float *)malloc(sizeof(float)*nT));
  assert(bE=(float *)malloc(sizeof(float)*nT));
  assert(E=(float *)malloc(sizeof(float)*nT));
  assert(D=(float *)malloc(sizeof(float)*nrow));
   
  /********************************************************/
  /***   compute E by permutation                       ***/
  /********************************************************/

  /* compute the original stat */
  Rprintf("\nE of the orginial data is: "); 
  for (i=0;i<nT;i++) {
    func_compute_stat=ptd->stat_array[i];
    (*func_compute_stat)(pdata, L, tmpT, &extras[i]);
    for (j=0;j<nrow;j++){
      pdr->T[j][i]=tmpT[j];
      if(func_max==max_abs) tmpT[j]=fabs(tmpT[j]);
      T[j][i]=tmpT[j];
    }
    E[i] = (*func_max)(tmpT, nrow);
    Rprintf("%5.3f  ", E[i]);
    if (weighted_dist) {
      wval=mad(tmpT, nrow);
      wval=1/(wval*wval);
      pdr->wval[i]=wval;
    }
    else pdr->wval[i]=1;
  }
  
   
  /* start permutation */
  creat_sampling(ncol, L, (*B));
  is_next=(*func_next_sample)(bL);
  b=0;
  while(is_next){
    for(i=0;i<nT;i++) {
      func_compute_stat = ptd->stat_array[i];
      (*func_compute_stat)(pdata, bL, tmpT, &extras[i]);
      bE[i]=(*func_max)(tmpT, nrow);
      if((func_max==max_high)&&(bE[i]>E[i]))
	E[i]=bE[i];
      else if ((func_max==max_low)&&(bE[i]<E[i]))
	E[i]=bE[i];
      else if ((func_max==max_abs)&&(fabs(bE[i])>fabs(E[i])))
	E[i]=bE[i];
    }
    b++;
    print_b(b,*B,"b=");
    is_next=(*func_next_sample)(bL);
  }

  Rprintf("\nAfter permutation , E is set at: "); 
  for (i=0;i<nT;i++) {
    Rprintf("%5.3f  ", E[i]);
    pdr->E[i]=E[i];
  }
  Rprintf("\n");
  
  /* compute distance for original stats */
  compute_euclid(T, nrow, nT, E, pdr->wval, pdr->D); 
  
  /* order the distance */
  order_index(pdr->D, pdr->R, nrow); /* R stores the index of sorted distance */
  
  free(bL);
  free(L);
  free(extras);
  free(tmpT);
  free(D);
  for (i=0;i<nrow;i++) free(T[i]);
  free(T);
  delete_sampling();
  
}     
    

void func_get_FDR(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B)
{
  int b=0, is_next; /* permutation index */
  int *L, *bL; /* L: class index label; bL: permutated class label */
  int i, j; /* looping index */
  float **bT, *tmpT; /* T: nrow by nT, stores original stats; bT: permutation stats */
  float *bD, *bMD; /* D: original distance vector; bD: permutation distance */
  float *fE, *fF;
  int ncol=pdata->ncol; /* number of samples */
  int nrow=pdata->nrow; /* number of genes */
  int nT=ptd->n_stat; /* number of statistics, t, fc, etc */
  int nsig=pdr->nsig;
  FUNC_COMPUTE_STAT func_compute_stat;
  FUNC_MAX func_max=ptd->func_max;
  FUNC_SAMPLE func_next_sample=ptd->func_next_sample;
  FUNC_COMPUTE_P func_compute_p=ptd->func_compute_p;
  float *extras;
  
  assert(extras=(float *)malloc(nT*sizeof(float)));
  memcpy(extras, ptd->extras, nT*sizeof(float));
  assert(L=(int *)malloc(sizeof(int)*ncol));
  memcpy(L, pdata->L, sizeof(int)*ncol);
  assert(bL=(int *)malloc(sizeof(int)*ncol));
  assert(tmpT=(float *)malloc(sizeof(float)*(nrow)));
  assert(bT=(float **)malloc(sizeof(float*)*nrow));
  for(i=0;i<nrow;i++)
    assert(bT[i]=(float *)malloc(sizeof(float)*nT));
  assert(fE=(float *)malloc(sizeof(float)*nT));
  assert(fF=(float *)malloc(sizeof(float)*nrow));
  assert(bD=(float *)malloc(sizeof(float)*nrow));
  assert(bMD=(float *)malloc(sizeof(float)*(nrow*(*B))));
    
  creat_sampling(ncol, L, (*B));
  is_next=(*func_next_sample)(bL);
  b=0;
  while(is_next){
     for (i=0;i<nT;i++) {
	func_compute_stat=ptd->stat_array[i];
	(*func_compute_stat)(pdata, bL, tmpT, &extras[i]);
	for (j=0;j<nrow;j++) {
	  bT[j][i]=tmpT[j];
	  if(func_max==max_abs) bT[j][i]=fabs(bT[j][i]);
      }
    }
    compute_euclid(bT, nrow, nT, pdr->E, pdr->wval, bD);
    /*qsort(bD, nrow, sizeof(bD[0]), distCompare);*/
    for(i=0;i<nrow;i++)  bMD[b*nrow+i]=bD[i];
        
    b++;
    print_b(b,*B,"b=");
    is_next=(*func_next_sample)(bL);
    }
 
  (*func_compute_p)(bMD, pdr->D, pdr->R, &nrow, B, &nsig, fF);
  for(i=0;i<nrow;i++) pdr->FDR[i]=(double)fF[i];
    
  free(tmpT);
  for (i=0;i<nrow;i++)
    free(bT[i]);
  free(bL);
  free(bD);
  free(bMD);
  free(fE);
  free(fF);
  delete_sampling();
}
      

void print_b(int b, int B, char *prompt) 
{
  static int p=0;
  if(b==0) p=0;
  if ((B<=100) || (b%(B/100)==0)) {
    Rprintf("%s%d\t", prompt,b);
    p++;
    if(p%10==0)
      Rprintf("\n");
  }
}


int type2test(char **options, TEST_DATA *td, int *nT, int *nL, float *extras)
{
  int i;
  assert(td->stat_array=(FUNC_COMPUTE_STAT *)malloc((*nT)*sizeof(FUNC_COMPUTE_STAT)));
  assert(td->extras=(float *)malloc((*nT)*sizeof(float)));
  
  if(*nL==1) Rprintf("\nOne-sample Statistics:\n");
  else if(*nL==2) Rprintf("\nTwo-sample Statistics:\n");
  else Rprintf("\nMulti-sample Statistics:\n");

  for (i=0;i<(*nT);i++) {
    (td->stat_array)[i]=type2stat(options[i], nL);
    (td->extras)[i]=extras[i];
    /*Rprintf("%5.3f ", td->extras[i]);*/
  }
  /*Rprintf("\n");*/
  
  td->n_stat=*nT;
  if(strcmp(options[*nT],"abs")==0){
    td->func_max=max_abs;
    td->func_cmp=cmp_abs;
  }
  else if(strcmp(options[*nT],"lower")==0){
    td->func_max=max_low;
    td->func_cmp=cmp_low;
  }
  else if(strcmp(options[*nT],"higher")==0){
    td->func_max=max_high;  
    td->func_cmp=cmp_high;
  }
  else 
    return 0;
  
  if(strcmp(options[*nT+1],"euclid")==0)
    td->weighted_dist=0;
  else if(strcmp(options[*nT+1],"weuclid")==0)
    td->weighted_dist=1;
  else
    return 0;

  if(strcmp(options[*nT+2],"fdr")==0)
    td->func_compute_p=calc_FDR;
  else if(strcmp(options[*nT+2],"adjp")==0)
    td->func_compute_p=calc_adjP;
  else
    return 0;

  if (*nL==1) td->func_next_sample=next_sample_1;
  else td->func_next_sample=next_sample;

  return 1;
}
  
  
FUNC_COMPUTE_STAT type2stat(char *ptest, int *nL)
{
  FUNC_COMPUTE_STAT test;
  if((strcmp(ptest,"t")==0) && (*nL==2)){
    Rprintf("t \t");
    test=compute_t2_stat;
   }
  else if((strcmp(ptest,"t")==0) && (*nL==1)){
    Rprintf("t \t");
    test=compute_t1_stat;
  }
  else if((strcmp(ptest,"fc")==0) && (*nL==1)){
    Rprintf("FC \t");
    test=compute_fc1_stat;
  }
  else if((strcmp(ptest,"fc")==0) && (*nL==2)){
    Rprintf("FC \t");
    test=compute_fc2_stat;
  }
  else if((strcmp(ptest,"fc")==0) && (*nL>2)) {
    Rprintf("FC \t");
    test=compute_fcm_stat;
  }
  else if((strcmp(ptest,"sam")==0) && (*nL==1)) {
    Rprintf("SAM \t");
    test=compute_sam1_stat_q;
  }
  else if((strcmp(ptest,"sam")==0) && (*nL==2)){
    Rprintf("SAM \t");
    test=compute_sam2_stat_q;
  }
  else if((strcmp(ptest,"sam")==0) && (*nL>2)) {
    Rprintf("SAM \t");
    test=compute_samm_stat_q;
  }
  else if((strcmp(ptest,"f")==0) && (*nL>=2)) {
    Rprintf("F \t");
    test=compute_f_stat;
  }
  else if((strcmp(ptest,"modt")==0) && ((*nL==2)||(*nL==1))){
    Rprintf("moderated t \t");
    test=compute_t_mod_stat;
   }
  else if((strcmp(ptest,"modf")==0) && (*nL>=2)){
    Rprintf("moderated F \t");
    test=compute_f_mod_stat;
   }
  else if((strcmp(ptest,"B")==0) && ((*nL==2)||(*nL==1))){
    Rprintf("B \n");
    test=compute_t_mod_B;
   }
  else error("invalid statistic parameter");
  return test;
}

/** bD -- distance matrix, ncol--no of permutation , nrow--no of genes **/
/** R -- order of D */
void calc_FDR(float *bD, float *D, int *R, int *pnrow, int *pncol, int *nsig, float *F)
{
  float **bMD, **count;
  int i, j, m;
    
  assert(bMD=(float **)malloc(sizeof(float *)*(*pnrow)));
  for(i=0;i<(*pnrow);i++)
    assert(bMD[i]=(float *)malloc(sizeof(float)*(*pncol)));
  assert(count=(float **)malloc(sizeof(float *)*(*nsig)));
  for(i=0;i<(*nsig);i++){
    assert(count[i]=(float *)malloc(sizeof(float)*(*pncol)));
    memset(count[i], 0, sizeof(float)*(*pncol));
  }

  for(i=0;i<(*pnrow);i++){
    for(j=0;j<(*pncol);j++){
      bMD[i][j]=bD[j*(*pnrow)+i];
    }
  }

  for(i=0;i<(*pncol);i++) {
    for(j=0;j<(*nsig);j++) { 
      int k=0;
      for(m=0;m<(*pnrow);m++) 
	if(bMD[m][i]<=D[R[j]]) k++;
      count[j][i]=(float)k;    
    }
  }
  
  for(i=0;i<(*nsig);i++){
    if(!R_FINITE(D[i])) F[i]=NA_REAL;
    else
      F[i]=median(count[i], (*pncol))/(float)(i+1);
      }

  /*for (i=0;i<(*nsig);i++) {
    if(!R_FINITE(D[i])) F[i]=NA_REAL;
    else {
      float m=0;
      for (j=0;j<(*pncol);j++)
	m+=count[i][j];
      F[i]=m/(float)(*pncol)/(float)(i+1);
    }
    }*/

  /* for(i=1;i<(*nsig);i++)
     if(F[i]<F[i-1]) F[i]=F[i-1];*/
  for(i=(*nsig-1);i>0;i--)
    if(F[i-1]>F[i]) F[i-1]=F[i];
  for(i=0;i<(*nsig);i++)
     if(F[i]>1.0) F[i]=1.0;
  for(i=(*nsig);i<(*pnrow);i++) F[i]=1.0;

 
  for(i=0;i<(*pnrow);i++) 
      free(bMD[i]);
  free(bMD);
  for(i=0;i<(*nsig);i++) 
   free(count[i]);
  free(count);
}

/*************************************************************************************/
/* pnrow -- numbe of genes (nG); pncol -- number of permutation (B)                 **/
/* bD -- dimension nGXB                                                             **/
/* R -- gene order                                                                  **/
/*************************************************************************************/
void calc_adjP(float *bD, float *D, int *R, int *pnrow, int *pncol, int *nsig, float *F)
{
  float **bMD, *Adj_P;
  int i, j, *count, *total;
    
  assert(bMD=(float **)malloc(sizeof(float *)*(*pnrow)));
  for(i=0;i<(*pnrow);i++)
    assert(bMD[i]=(float *)malloc(sizeof(float)*(*pncol)));
  assert(count=(int *)malloc(sizeof(int)*(*pnrow)));
  assert(total=(int *)malloc(sizeof(int)*(*pnrow)));
  assert(Adj_P=(float *)malloc(sizeof(float)*(*pnrow)));
  memset(count, 0, sizeof(int)*(*pnrow));
  memset(total, 0, sizeof(int)*(*pnrow));

  for(i=0;i<(*pnrow);i++){
    for(j=0;j<(*pncol);j++){
      bMD[i][j]=bD[j*(*pnrow)+i];
    }
  }
   
  for(i=0;i<(*pncol);i++) {
     float qT=bMD[R[*pnrow-1]][i];/*intitalize the qT*/
     if(qT<=D[R[*pnrow-1]]) count[*pnrow-1]+=1;
     if(R_FINITE(qT)) total[*pnrow-1]+=1;
     for(j=*pnrow-2;j>=0;j--){ /*looping the row reversely*/
       if(!R_FINITE(D[R[j]])) continue;
       if(bMD[R[j]][i]<qT) qT=bMD[R[j]][i];
       if((R_FINITE(bMD[R[j]][i]))&&(!R_FINITE(qT)))
	  qT=bMD[R[j]][i];
       if(qT<=D[R[j]]) count[j]+=1;
       if(R_FINITE(qT)) total[j]+=1;
     }
  }

  for(i=0;i<(*pnrow);i++){ 
    if(total[i]==0) 
      Adj_P[i]=NA_REAL;
    else Adj_P[i]=count[i]*1.0/total[i];
  }
 
  /*enforce the montonicity*/
  for(i=1;i<(*pnrow);i++)
    if(Adj_P[i]<Adj_P[i-1])
      Adj_P[i]=Adj_P[i-1];

  for(i=0;i<(*pnrow);i++)
    F[i]=Adj_P[i];

  for(i=0;i<(*pnrow);i++) free(bMD[i]);
  free(bMD);
  free(count);
  free(total);
  free(Adj_P);

}


void extract_deds_res(DEDS_RES *pdr, double *E, int *R, double *FDR, double *T)
{
  int nr=pdr->nrow;
  int nT=pdr->nT;
  int i,j;

  for(i=0;i<nT;i++)
    E[i]=(double)pdr->E[i];
  for(i=0;i<nr;i++) {
    R[i]=pdr->R[i];
    FDR[i]=pdr->FDR[i];
  }
  for(i=0;i<nr;i++)
    for(j=0;j<nT;j++)
      T[nr*j+i]=(double)pdr->T[i][j];
 
}
    

