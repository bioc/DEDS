#include "utilities.h"


void get_stat(double *d, int *pnrow, int *pncol, int *L, float *T, char **options, float *extras, int *nL)
{
  GENE_DATA data;
  FUNC_COMPUTE_STAT test;
 
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  test=type2stat(options[0], nL);
  (*test)(&data, data.L, T, &extras[0]);
   free_gene_data(&data);
}

void get_unadjp(double *d, int *pnrow, int *pncol, int *L, float *T, float *P, char **options, float *extras, int *nL, int *B)
{
  GENE_DATA data;
  TEST_DATA td;
  FUNC_COMPUTE_STAT func_test;
  FUNC_MAX func_max;
  FUNC_SAMPLE func_next_sample;
  int nT=1; /* handling one test */
  float *TB;
  int b, *bL, is_next;
  int i;
  int *count, *total;
  
  assert(TB=(float *)malloc(sizeof(float)*(*pnrow)));
  assert(bL=(int *)malloc(sizeof(int)*(*pncol)));
  assert(count=(int *)malloc(sizeof(int)*(*pnrow)));
  memset(count, 0, sizeof(int)*(*pnrow));
  assert(total=(int *)malloc(sizeof(int)*(*pnrow)));
  memset(total, 0, sizeof(int)*(*pnrow));
  

  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  if(type2test(options, &td, &nT, nL, extras)==0) return;
  func_max=td.func_max;
  func_next_sample=td.func_next_sample;
  func_test=td.stat_array[0];
 
  (*func_test)(&data, data.L, T, &extras[0]);
  
  creat_sampling(*pncol, L, (*B));
  is_next=func_next_sample(bL);
  b=0;
 
  while(is_next){
    (*func_test)(&data, bL, TB, nL);
    for(i=0;i<(*pnrow);i++){
      if(!R_FINITE(T[i])) continue;
      if(!R_FINITE(TB[i])) continue;
      if((func_max==max_high)&&(TB[i]>=T[i]))
	  count[i]++;
      else if ((func_max==max_low)&&(TB[i]<=T[i]))
	  count[i]++;
      else if ((func_max==max_abs)&&(fabs(TB[i])>=fabs(T[i])))
	  count[i]++;
      total[i]++;
    }
    b++;
    print_b(b, *B, "b=");
    is_next=func_next_sample(bL);
  }
  
  for(i=0;i<(*pnrow);i++) {
    if(total[i]==0) P[i]=NA_REAL;
    else P[i]=(double)(count[i]*1.0)/(total[i]*1.0);
    }
  
  free(count);
  free(total);
  free(TB);
  free(bL);
  free_gene_data(&data);
  delete_sampling();
}

/********************************************************************************/
/*  modified from multtest package function adj-by_t                                                          */
/********************************************************************************/
/*  
  L: is the labelling of each experiment
  T: is that test statistics
  P: unadjtesed P-values
  Adj_P:ajusted p-values by using the max|T|
  To use the function first_sample, and next_sample, they're needed to write into a separate 
  file, where it provides the create_sampling to do initialization such as allocate the space
  (before use the sampling) and delete_sampling after we've done the sampling in the main().

  int first_sample(int *L)
       get the first sample of the labelling.
  int next_sample(int* L)
       get the next sample, if it's done all the sampling, then it returns 0, 
       otherwise it returns 1.


  input: d, L, B,options
  output: R,T,P,Adj_P
 */

void get_adjp(double *d, int *pnrow, int *pncol, int *L, float *T, float *P, float *Adj_P, int *R, char **options, float *extras, int *nL, int *B)
{
  GENE_DATA data;
  TEST_DATA td;
  FUNC_COMPUTE_STAT func_test;
  FUNC_CMP func_cmp;
  FUNC_SAMPLE func_next_sample;
  int nT=1; /* handling one test */
  float *BT, *count1, *count2, QT; /*QT successive maxima*/
  int b, *bL, is_next, *total1, *total2;
  int i;

  assert(BT=(float *)malloc(sizeof(float)*(*pnrow)));
  assert(bL=(int *)malloc(sizeof(int)*(*pncol)));
  assert(count1=(float*)malloc(sizeof(float)*(*pnrow)));
  memset(count1,0,sizeof(float)*(*pnrow)); 
  assert(total1=(int*)malloc(sizeof(int)*(*pnrow)));
  memset(total1,0,sizeof(int)*(*pnrow));
  assert(count2=(float*)malloc(sizeof(float)*(*pnrow)));
  memset(count2,0,sizeof(float)*(*pnrow)); 
  assert(total2=(int*)malloc(sizeof(int)*(*pnrow)));
  memset(total2,0,sizeof(int)*(*pnrow));

   /*comuter the original t-statfirst*/
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  if(type2test(options, &td, &nT, nL, extras)==0) return;
  func_cmp=td.func_cmp;
  func_next_sample=td.func_next_sample;
  func_test=td.stat_array[0];

  (*func_test)(&data, data.L, T, &extras[0]);

  /*sort the T*/  
  order_data(T,R,*pnrow,func_cmp);
  sort_gene_data(&data,R);
  sort_vector(T,R,*pnrow);
  
  /*iteration for permutaion*/
  creat_sampling(*pncol, L, (*B));
  is_next=func_next_sample(bL);
  b=0;

  /*changed to the orignal stat, which is monotone of t and centered*/
  while(is_next){
     (*func_test)(&data, bL, BT, nL);
    /*deal with unajdused value first*/
    for(i=0;i<(*pnrow);i++){
      if(T[i]==NA_FLOAT) continue;
      if(BT[i]!=NA_FLOAT){
	if((func_cmp==cmp_high)&&(BT[i]+EPSILON>=T[i])) count2[i]++;
	if((func_cmp==cmp_low)&&(BT[i]<=T[i]+EPSILON)) count2[i]++;
	if((func_cmp==cmp_abs)&&(fabs(BT[i])>=fabs(T[i])-EPSILON)) count2[i]++;
	total2[i]++;
      }
    }

    /*deal with adjusted values*/
    QT=NA_FLOAT;/*intitalize the QT*/
    for(i=(*pnrow-1);i>=0;i--){ /*looping the row reversely*/
      if(T[i]==NA_FLOAT) continue;
        /* right now I only implements the 3 cases, which are pretty common*/
      if(func_cmp==cmp_high){
	if((BT[i]!=NA_FLOAT)&&(QT!=NA_FLOAT)&&(BT[i]>QT))
	  QT=BT[i];
	if((BT[i]!=NA_FLOAT)&&(QT==NA_FLOAT))
	  QT=BT[i];
	if((QT!=NA_FLOAT)&&(QT>=T[i]-EPSILON)) count1[i]+=1;
      }else if(func_cmp==cmp_low){
	if((BT[i]!=NA_FLOAT)&&(QT!=NA_FLOAT)&&(BT[i]<QT))
	  QT=BT[i];
	if((BT[i]!=NA_FLOAT)&&(QT==NA_FLOAT))
	  QT=BT[i];
	if((QT!=NA_FLOAT)&&(QT<=T[i]+EPSILON)) count1[i]+=1;
      }else if (func_cmp==cmp_abs) {
	if((BT[i]!=NA_FLOAT)&&(QT!=NA_FLOAT)&&(fabs(BT[i])>QT))
	  QT=fabs(BT[i]);
	if((BT[i]!=NA_FLOAT)&&(QT==NA_FLOAT))
	  QT=fabs(BT[i]);
	if((QT!=NA_FLOAT)&&(QT>=fabs(T[i])-EPSILON)) count1[i]+=1;
      }	      
      if(QT!=NA_FLOAT) total1[i]++;
    }
    b++;
    print_b(b,*B,"b=");
    is_next=func_next_sample(bL);
  }

  /*summarize the results*/
  /*unadjusted one*/
  for(i=0;i<(*pnrow);i++){ 
    if(total2[i]==0) 
      P[i]=NA_FLOAT;
    else P[i]=count2[i]*1.0/total2[i];
  }
  /*adjused one*/
  for(i=0;i<(*pnrow);i++){ 
    if(total1[i]==0) 
      Adj_P[i]=NA_FLOAT;
    else Adj_P[i]=count1[i]*1.0/total1[i];
  }
  /*enforce the montonicity*/
  for(i=1;i<(*pnrow);i++)
    if(Adj_P[i]<Adj_P[i-1])
      Adj_P[i]=Adj_P[i-1];

  /*free the spaces*/
  free(BT);
  free(count1);
  free(total1);
  free(count2);
  free(total2);
  free(bL);
}      


void get_fdr(double *d, int *pnrow, int *pncol, int *L, float *T, float *P, float *FDR, int *R, char **options, float *extras, int *nL, int *B)
{
  GENE_DATA data;
  TEST_DATA td;
  FUNC_COMPUTE_STAT func_test;
  FUNC_CMP func_cmp;
  FUNC_SAMPLE func_next_sample;
  int nT=1; /* handling one test */
  float **BT, *junk, pi0, q1, q2; 
  int b, *bL, is_next, *R1, *total2, total1=0, *count1, *count2;
  int i,j;
  int total;

  assert(BT=(float **)malloc(sizeof(float *)*(*B)));
  for (i=0;i<(*B);i++) {
    assert(BT[i]=(float *)malloc(sizeof(float)*(*pnrow)));
  }
  assert(bL=(int *)malloc(sizeof(int)*(*pncol)));
  assert(count1=(int*)malloc(sizeof(int)*(*pnrow)));
  memset(count1,0,sizeof(int)*(*pnrow)); 
  assert(total2=(int*)malloc(sizeof(int)*(*pnrow)));
  memset(total2,0,sizeof(int)*(*pnrow));
  assert(count2=(int*)malloc(sizeof(int)*(*pnrow)));
  memset(count2,0,sizeof(int)*(*pnrow)); 
  
   /*compute the original t-statfirst*/
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  if(type2test(options, &td, &nT, nL, extras)==0) return;
  func_cmp=td.func_cmp;
  func_next_sample=td.func_next_sample;
  func_test=td.stat_array[0];

  (*func_test)(&data, data.L, T, &extras[0]);

  /*sort the T*/  
  order_data(T,R,*pnrow,func_cmp);
  sort_gene_data(&data,R);
  sort_vector(T,R,*pnrow);
  
  /*iteration for permutaion*/
  creat_sampling(*pncol, L, (*B));
  is_next=func_next_sample(bL);
  b=0;
  total=(*pnrow)*(*B);

  /*changed to the orignal stat, which is monotone of t and centered*/
  while(is_next){
    (*func_test)(&data, bL, BT[b], nL);
    /*deal with unajdused value first*/
      for(i=0;i<(*pnrow);i++){
      if(T[i]==NA_FLOAT) continue;
      if(BT[b][i]!=NA_FLOAT){
	if((func_cmp==cmp_high)&&(BT[b][i]+EPSILON>=T[i])) count2[i]++;
	if((func_cmp==cmp_low)&&(BT[b][i]<=T[i]+EPSILON)) count2[i]++;
	if((func_cmp==cmp_abs)&&(fabs(BT[b][i])>=fabs(T[i])-EPSILON)) count2[i]++;
	total2[i]++;
      }
    }
    b++;
    print_b(b,*B,"b=");
    is_next=func_next_sample(bL);
  }

  assert(junk=(float*)malloc(sizeof(float)*total));
  assert(R1=(int*)malloc(sizeof(int)*total));
  for (i=0;i<(*B);i++){
    for (j=0;j<(*pnrow);j++){
      junk[i*(*pnrow)+j]=BT[i][j];
    }
  }
   
   
  order_data(junk,R1,total,func_cmp);
  sort_vector(junk,R1,total);
	 
  /*deal with fdr*/
  for(i=0;i<(*pnrow);i++){
    if(T[i]==NA_FLOAT) continue;
    for(j=0;j<total;j++) {
      if(junk[j]==NA_FLOAT) continue;
      if((func_cmp==cmp_high) && (junk[j]<T[i])) {
	count1[i]=j;
	break;
      }
      if((func_cmp==cmp_low) && (junk[j]>T[i])) {
	count1[i]=j;
	break;
      }
      if((func_cmp==cmp_abs) && (fabs(junk[j])<fabs(T[i]))) {
	count1[i]=j;
	break;
      }
     }
    if(count1[i]>(i+1)*(*B)) {
      for(j=(i+1);j<(*pnrow);j++) count1[j]=(j+1)*(*B);
      break;
    }
  }
  
  q1=T[(int)floor(0.25*(*pnrow))];
  q2=T[(int)floor(0.75*(*pnrow))];
  for (i=0;i<total;i++) {
    if((func_cmp==cmp_high) && (junk[i]<q1) && (junk[i]>q2)) total1++;
    if((func_cmp==cmp_low) && (junk[i]>q1) && (junk[i]<q2)) total1++;
    if((func_cmp==cmp_abs) && (fabs(junk[i])<fabs(q1)) && (fabs(junk[i])>fabs(q2))) total1++;
  }
  pi0=total1*1.0/(*B)/(0.5*(*pnrow));
  if(pi0>1) pi0=1;
  Rprintf("\npi0=%5.3f\n", pi0);
    
  /*summarize the results*/
  /*unadjusted one*/
  for(i=0;i<(*pnrow);i++){ 
    if(total2[i]==0) 
      P[i]=NA_FLOAT;
    else P[i]=count2[i]*1.0/total2[i];
  }
  for(i=0;i<(*pnrow);i++)
    FDR[i]=pi0*(count1[i]*1.0)/(*B)/(i+1);
  
  for(i=1;i<(*pnrow);i++){
    if(FDR[i]<FDR[i-1])
      FDR[i]=FDR[i-1];
  }

  for(i=0;i<(*pnrow);i++) 
    if(FDR[i]>1) FDR[i]=1;

  /*free the spaces*/
  free(BT);
  free(count1);
  free(total2);
  free(count2);
  free(junk);
  free(bL);
}      


/*******************************************************************/
/**                    t stats                                    **/
/*******************************************************************/

void compute_t2_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i;
  for(i = 0; i < pdata->nrow; i++)
    T[i]= t2_stat(pdata->d[i], L, pdata->ncol, extra);
  
}
float t2_stat(const float *Y, const int *L, const int n, const void *extra)
{
  float mean_na[2]={0,0}, ss_na[2]={0,0}, c0, c1, num, denum;
  int i, count[2]={0,0}, class;
     
  /* calculating mean and number of non-na values */
  for (i = 0; i < n; i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class = L[i];
    mean_na[class] += Y[i];
    count[class]++;
  }

  mean_na[0] /= (count[0] * 1.0);
  mean_na[1] /= (count[1] * 1.0);

  /* computing variance */
  for (i = 0; i < n; i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class = L[i];
    ss_na[class] += ((Y[i] - mean_na[class])) * ((Y[i]-mean_na[class]));
  }
  
  if(ss_na[0] + ss_na[1] == 0) return NA_REAL;
  c0 = (count[0] * (count[0] - 1));
  c1 = (count[1] * (count[1] - 1));
  num = mean_na[0] - mean_na[1];
  denum = sqrt(ss_na[0]/c0 + ss_na[1]/c1);
  return(num/denum);
}

void compute_t1_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i, k=1;
  for (i=0;i<pdata->ncol;i++) {
    if (L[i]==-1) {
      k=-1;
      break;
    }
  }
  for(i = 0; i < pdata->nrow; i++)
    T[i]= t1_stat(pdata->d[i], L, pdata->ncol, &k); 
}

float t1_stat(const float *Y, const int *L, const int n, const void *extra)
{
  float mean_na=0, ss_na=0;
  int i, count=0;
  int k=*(int *)extra;
  
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    if (k==1) mean_na+=Y[i];
    else mean_na+=L[i]*Y[i];
    count++;
  }
  mean_na/=(count*1.0);

  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    if (k==1) ss_na+=(Y[i]-mean_na)*(Y[i]-mean_na);
    else ss_na+=(L[i]*Y[i]-mean_na)*(L[i]*Y[i]-mean_na);
  }
  ss_na/=(float)(count-1);
  ss_na=sqrt(ss_na/(float)count);
  
  if (ss_na==0) return NA_REAL;
  else return(mean_na/ss_na);
}


/*******************************************************************/
/**                    FC                                         **/
/*******************************************************************/

void compute_fc2_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i;
  for(i = 0; i < pdata->nrow; i++)
    T[i]= fc2_stat(pdata->d[i], L, pdata->ncol, extra); 
}

float fc2_stat(const float *Y, const int *L, const int n, const void *extra)
{
  float mean_na[2]={0,0};
  int i, count[2]={0,0}, class;
 
  /* calculating mean and number of non-na values */
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class=L[i];
    mean_na[class]+=Y[i];
    count[class]++;
  }

  if ((count[0]==0) || (count[1]==0)) return(NA_REAL);
  mean_na[0]/=(count[0]*1.0);
  mean_na[1]/=(count[1]*1.0);
  return(mean_na[0]-mean_na[1]);
}

 
void compute_fc1_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i, k=1;
  for (i=0;i<pdata->ncol;i++) {
    if (L[i]==-1) {
      k=-1;
      break;
    }
  }
  for(i = 0; i < pdata->nrow; i++)
    T[i]= fc1_stat(pdata->d[i], L, pdata->ncol, &k); 
}

float fc1_stat(const float *Y, const int *L, const int n, const void *extra)
{
  float mean_na=0;
  int i,count=0;
  int k=*(int *)extra;
  
  
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    if (k==1) mean_na+=Y[i];
    else mean_na+=L[i]*Y[i];
    count++;
  }
  if (count==0) return(NA_REAL); 
  mean_na/=(float)count;
  return(mean_na);
}

void compute_fcm_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i;
  for(i = 0; i < pdata->nrow; i++)
    T[i]= fcm_stat(pdata->d[i], L, pdata->ncol, &(pdata->nL)); 
}

float fcm_stat(const float *Y, const int *L, const int n,  const void *extra)
{
  float *mean_na, max_m, min_m;
  int i, class, *count;
  int nL=*(int *)extra;
 
  assert(mean_na=(float *)malloc(nL*sizeof(float)));
  memset(mean_na, 0, nL*sizeof(float));
  assert(count=(int *)malloc(nL*sizeof(int)));
  memset(count, 0, nL*sizeof(int));
    
  /* calculating mean and number of non-na values */
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class=L[i];
    mean_na[class]+=Y[i];
    count[class]++;
  }
  
  for (i=0;i<nL;i++){
    if (count[i]==0) return(NA_REAL);
    mean_na[i]/=(count[i]*1.0);
  }

  max_m=max_high(mean_na, nL);
  min_m=max_low(mean_na, nL);
  return(max_m-min_m);
}


/*******************************************************************/
/**                    SAM stats                                  **/
/*******************************************************************/
void compute_sam2_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  float *num, *denum, s0, c0, c1;
  int i, j, class;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  float q=*(float *)extra;
  
  assert(num=(float *)malloc(sizeof(float)*nrow));
  assert(denum=(float *)malloc(sizeof(float)*nrow));
  
  /* calculating mean and number of non-na values */
  for (i=0;i<nrow;i++) {
    float mean_na[2]={0,0}, ss_na[2]={0,0};
    int count[2]={0,0};
    
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j]))  continue;
      class = L[j];
      mean_na[class]+=pdata->d[i][j] ;
      count[class]++;
    }
    mean_na[0]/=(count[0]*1.0);
    mean_na[1]/=(count[1]*1.0);
        
    /* computing variance */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      class=L[j];
      ss_na[class]+=(pdata->d[i][j]-mean_na[class])*(pdata->d[i][j]-mean_na[class]);
    }
    
    if((ss_na[0]==0)||(ss_na[1]==0)) denum[i]=NA_REAL;
    else {
      c0=1/(1.0*count[0])+1/(1.0*count[1]);
      c1=(float)(count[0]+count[1] - 2);
      num[i]=mean_na[0]-mean_na[1];
      denum[i]=sqrt((ss_na[0] + ss_na[1])*c0/c1);
     }
  }
  quantile(denum, nrow, &q, 1, &s0);
  
  for (i=0;i<nrow;i++) {
    if (denum[i]==NA_REAL) T[i]=NA_REAL;
    else T[i]=num[i]/(s0+denum[i]);
  }

  free(denum);
  free(num);
}

void compute_sam1_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  float *num, *denum, s0;
  int i, j;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  int k=1;  
  float q=*(float *)extra;
  
  assert(num=(float *)malloc(sizeof(float)*nrow));
  assert(denum=(float *)malloc(sizeof(float)*nrow));

  for (i=0;i<pdata->ncol;i++) {
    if (L[i]==-1) {
      k=-1;
      break;
    }
  }
  
  /* calculating mean and number of non-na values */
  for (i=0;i<nrow;i++) {
    float mean_na=0, ss_na=0;
    int count=0;
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j]))  continue;
      if (k==1) mean_na+=pdata->d[i][j] ;
      else mean_na+=L[j]*pdata->d[i][j] ;
      count++;
    }
    mean_na/=(count*1.0);
       
    /* computing variance */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      if (k==1) ss_na+=(pdata->d[i][j]-mean_na)*(pdata->d[i][j]-mean_na);
      else ss_na+=(L[j]*pdata->d[i][j]-mean_na)*(L[j]*pdata->d[i][j]-mean_na);
    }
        
    if(count==0) denum[i]=NA_REAL;
    else {
      num[i]=mean_na;
      denum[i]=sqrt(ss_na/(count*(count-1)*1.0));
    }
  }
  quantile(denum, nrow, &q, 1, &s0);
  for (i=0;i<nrow;i++) {
    if (denum[i]==NA_REAL) T[i]=NA_REAL;
    else T[i]=num[i]/(s0+denum[i]);
  }

  free(denum);
  free(num);
}

void compute_samm_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  float *num, *denum, s0, dev;
  int i, j,  class;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  int nL=pdata->nL;
  float q=*(float *)extra;
  
  assert(num=(float *)malloc(sizeof(float)*nrow));
  assert(denum=(float *)malloc(sizeof(float)*nrow));

  /* calculating mean and number of non-na values */
  for (i=0;i<nrow;i++) {
    float meani[nL], ssi[nL], mean=0, N1=1.0, N2=0.0;
    int counti[nL], N=0;
    float wss=0, bss=0;
    memset(meani, 0, sizeof(float)*nL);
    memset(ssi, 0, sizeof(float)*nL);
    memset(counti, 0, sizeof(int)*nL);
    
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      class=L[j];
      meani[class]+=pdata->d[i][j];
      counti[class]++;
      N++;
      mean+=pdata->d[i][j];
    }

    mean/=(N*1.0);
    for (j=0;j<nL;j++)  meani[j]/=(counti[j]*1.0);
       
    for (j=0;j<ncol;j++) {
      if(!R_FINITE(pdata->d[i][j])) continue;
      class=L[j];
      dev=pdata->d[i][j]-meani[class];
      ssi[class]+=dev*dev;
    }
    
    for (j=0;j<nL;j++) {
      wss+=ssi[j];
      dev=meani[j]-mean;
      bss+=dev*dev*counti[j];
      
      N1*=counti[j];
      N2+=(1/(1.0*counti[j]));
    }
    
    num[i]=(N/N1)*bss/(nL-1.0);
    denum[i]=N2*wss/(N-nL-0.0);
  }
  
  quantile(denum, nrow, &q, 1, &s0);
  for (i=0;i<nrow;i++) T[i]=num[i]/(denum[i]+s0);

  free(num);
  free(denum);
}  
  
/*******************************************************************/
/**                    F stats                                    **/
/*******************************************************************/

void compute_f_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  int i;
  for(i = 0; i < pdata->nrow; i++)
    T[i]= f_stat(pdata->d[i], L, pdata->ncol, &(pdata->nL)); 
}
      
float f_stat(const float *Y, const int *L, const int n, const void *extra)
{
  int nL=*(int *)extra, valid=1;
  float mean=0, meani[nL], ssi[nL], num, denum, bss=0, wss=0, dev;
  int i, counti[2], class, N=0;
    
  memset(meani, 0, sizeof(float)*nL);
  memset(ssi, 0, sizeof(float)*nL);
  memset(counti, 0, sizeof(int)*nL);
  
  /* calculating mean and number of non-na values */
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class=L[i];
    meani[class]+=Y[i];
    counti[class]++;
    mean+=Y[i];
    N++;
  }

  mean/=(N*1.0);
 
  for (i=0;i<nL;i++)   {
    if ((counti[i]==0)||(counti[i]==1)) {
      valid=0;
      break;
    }
    else meani[i] /= (counti[i] * 1.0);
  }
  
  /* computing variance */
  for (i=0;i<n;i++) {
    if (!R_FINITE(Y[i]))
      continue;
    class = L[i];
    ssi[class]+=((Y[i]-meani[class]))*((Y[i]-meani[class]));
  }
  
  for (i=0;i<nL;i++) {
    if (ssi[i]==0) {
      valid=0;
      break;
    }
    wss+=ssi[i];
    dev=meani[i]-mean;
    bss+=dev*dev*counti[i];
  }
  
  if(valid==0) return NA_REAL;
  else {
    num=bss/(nL-1.0);
    denum=wss/(N-nL-0.0);
    return(num/denum);
  }  
}

/*******************************************************************/
/**                  moderated t                                  **/
/*******************************************************************/
/*void F77_NAME(ch2inv)(double *x, int *ldx, int *n, double *v, int *info);
  void F77_NAME(chol)(double *a, int *lda, int *n, double *v, int *info);*/

void inv_chol(double *x, int *nr, int *nc, double *v, int *info)
{
  double *dx;
  int n=*nr;
  
  if((*nr)!=(*nc)) error("non-square matrix");
  assert(dx=(double *)malloc(sizeof(double)*n*n));
  memset(dx, 0, sizeof(double)*n*n);

  /* do choleski decomposition */
  F77_CALL(chol)(x, &n, &n, dx, info);
  /* do inverse of choleski decomposition */
  F77_CALL(ch2inv)(dx, &n, &n, v, info);
 
}

void t2_mod_stat_func(GENE_DATA *pdata, int *L, TMOD_DATA *ptmod)
{
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  int nL=pdata->nL;
  int i,j;
  
  for (i=0;i<nrow;i++) {
    float mean[2]={0,0}, ss[2]={0,0};
    int count[2]={0,0}, class;
        
    /* calculating mean and number of non-na values */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j]))  continue;
      class = L[j];
      mean[class]+=pdata->d[i][j] ;
      count[class]++;
    }
    mean[0]/=(count[0]*1.0);
    mean[1]/=(count[1]*1.0);
    
    ptmod->stdev_unscale[i]=sqrt(1.0/count[0]+1.0/count[1]);
 
    /* computing variance */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      class=L[j];
      ss[class]+=(pdata->d[i][j]-mean[class])*(pdata->d[i][j]-mean[class]);
    }
    
    ptmod->df_resid[i]=count[0]+count[1]-nL;
    ptmod->sigma2[i]=(ss[0]+ss[1])/ptmod->df_resid[i];
    ptmod->mean[i]=mean[0]-mean[1];
  }
}
   

void t1_mod_stat_func(GENE_DATA *pdata, int *L, TMOD_DATA *ptmod)
{
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  int nL=pdata->nL;
  int i,j;
 
  for (i=0;i<nrow;i++) {
    float m=0, ss=0;
    int count=0;
        
    /* calculating mean and number of non-na values */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j]))  continue;
      m+=pdata->d[i][j] ;
      count++;
    }
    m/=(count*1.0);
   
    /* computing variance */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      ss+=(pdata->d[i][j]-m)*(pdata->d[i][j]-m);
    }
    ptmod->df_resid[i]=count-nL;
    ptmod->sigma2[i]=ss/ptmod->df_resid[i];
    ptmod->mean[i]=m;
    ptmod->stdev_unscale[i]=sqrt(1.0/count);
  }
}


void compute_t_mod_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  TMOD_DATA tmod;
  float s2_prior=0, df_prior=0, *s2_post;
  int i;
  int nrow=pdata->nrow;
  int nL=pdata->nL;
   
  create_tmod_data(&nrow, &tmod);
  assert(s2_post=(float *)malloc(sizeof(float)*nrow));
  if(nL==1) t1_mod_stat_func(pdata, L, &tmod);
  else t2_mod_stat_func(pdata, L, &tmod);
    
  fitFDist(tmod.sigma2, tmod.df_resid, nrow, &s2_prior, &df_prior);
  
  for(i=0;i<nrow;i++) {
    /* calculate s2 posterior */
    if(!R_FINITE(df_prior)) s2_post[i]=s2_prior;
    else if(tmod.df_resid[i]==0) s2_post[i]=(df_prior*s2_prior)/df_prior;
    else s2_post[i]=((float)tmod.df_resid[i]*tmod.sigma2[i]+df_prior*s2_prior)/((float)tmod.df_resid[i]+df_prior);
    /* calculate moderated t */
    if(!R_FINITE(tmod.stdev_unscale[i])) T[i]=NA_REAL;
    else T[i]=tmod.mean[i]/(tmod.stdev_unscale[i]*sqrt(s2_post[i]));
  }
     
  free(s2_post);
  free_tmod_data(&tmod);
}


void get_t_mod_stat(double *d, int *pnrow, int *pncol, int *L, float *T, int *nL)
{
  GENE_DATA data;
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  compute_t_mod_stat(&data, L, T, nL);
  free_gene_data(&data);
 
}

void get_f_mod_stat(double *d, int *pnrow, int *pncol, int *L, float *T, int *nL)
{
  GENE_DATA data;
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  if(*nL<2) error("only one class, F tests can't be carried out, try t tests...");
  else compute_f_mod_stat(&data, L, T, nL);
  free_gene_data(&data);
}
 
void compute_f_mod_stat(GENE_DATA *pdata, int *L, float *T, const void *extra)
{
  float *sigma2, s2_prior=0, df_prior=0, *s2_post;
  float *bss, *wss;
  int i, j;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;
  int nL=pdata->nL;
  int *df_resid;
  
  assert(sigma2=(float *)malloc(sizeof(float)*nrow));
  assert(df_resid=(int *)malloc(sizeof(int)*nrow));
  assert(s2_post=(float *)malloc(sizeof(float)*nrow));
  assert(bss=(float *)malloc(sizeof(float)*nrow));
  memset(bss, 0, sizeof(float)*nrow);
  assert(wss=(float *)malloc(sizeof(float)*nrow));
  memset(wss, 0, sizeof(float)*nrow);
  
  for (i=0;i<nrow;i++) {
    float mean[nL], ss[nL], meanT=0;
    int count[nL], class, total=0;
    
    memset(mean, 0, sizeof(float)*nL);
    memset(ss, 0, sizeof(float)*nL);
    memset(count, 0, sizeof(int)*nL);
        
    /* calculating mean and number of non-na values */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j]))  continue;
      class = L[j];
      mean[class]+=pdata->d[i][j] ;
      count[class]++;
      meanT+=pdata->d[i][j];
      total++;
    }
    for (j=0;j<nL;j++) 
      mean[j]/=(count[j]*1.0);
    meanT=meanT/total;
    
    /* computing variance */
    for (j=0;j<ncol;j++) {
      if (!R_FINITE(pdata->d[i][j])) continue;
      class=L[j];
      ss[class]+=(pdata->d[i][j]-mean[class])*(pdata->d[i][j]-mean[class]);
    }
    
    for(j=0;j<nL;j++) {
      wss[i]+=ss[j];
      bss[i]+=(mean[j]-meanT)*(mean[j]-meanT)*count[j];
    }
    
    df_resid[i]=total-nL;
    sigma2[i]=wss[i]/df_resid[i];
  }
    
  fitFDist(sigma2, df_resid, nrow, &s2_prior, &df_prior);
  
  for(i=0;i<nrow;i++) {
    /* calculate s2 posterior */
    if(!R_FINITE(df_prior)) s2_post[i]=s2_prior;
    else if(df_resid[i]==0) s2_post[i]=(df_prior*s2_prior)/df_prior;
    else s2_post[i]=((float)df_resid[i]*sigma2[i]+df_prior*s2_prior)/((float)df_resid[i]+df_prior);
    /* calculate moderated f */
    if(!R_FINITE(s2_post[i])) T[i]=NA_REAL;
    else T[i]=bss[i]/(nL-1.0)/s2_post[i];
  }

  free(sigma2);
  free(df_resid);
  free(s2_post);
  free(bss);
  free(wss);
}



void fitFDist(float *sigma2, int *df1, int n, float *scale, float *df2)
{
  float *z, *e, emean=0, evar=0;
  int i, count=0;
  assert(z=(float *)malloc(sizeof(float)*n));
  assert(e=(float *)malloc(sizeof(float)*n));
  
  for(i=0;i<n;i++) {
     if((R_FINITE(sigma2[i]))&&(sigma2[i]>EPSILON)) {
       z[i]=log(sigma2[i]);
       e[i]=z[i]-digamma(df1[i]*1.0/2)+log(df1[i]*1.0/2);
       emean+=e[i];
       count++;
     }
  }
  emean=emean/(float)count;
  
  for(i=0;i<n;i++) {
    if((R_FINITE(sigma2[i]))&&(sigma2[i]>EPSILON))
      evar+=(float)count/(float)(count-1)*(e[i]-emean)*(e[i]-emean)-trigamma(df1[i]*1.0/2);
  }
  evar=evar/(float)count;
 
  if(evar>0) {
    *df2=2*trigammaInverse(evar);
    *scale=exp(emean+digamma(*df2/2)-log(*df2/2));
  }
  else {
    *df2=FLT_MAX;
    *scale=exp(emean);
  }
}

float trigammaInverse(float x)
{
  int iter;
  float tri, dif, y;

  if(x>1e+7) y=1/sqrt(x);
  else if(x<1e-6) y=1/x;
  else {
    y=0.5+1/x;
    for(iter=0;iter<50;iter++) {
      tri=trigamma(y);
      dif=tri*(1 - tri/(x))/tetragamma(y);
      y=y+dif;
      if(-dif/(y)<1e-08) break;
    }
  }
  return(y);
}

void testtmixture(float *t, int *n, float *std, float *df, float *proportion, float *c0lim, float *ret)
{
  ret[0]=tmixture(t, *n, std, df, *proportion, *c0lim);
}


float tmixture(float *t, int n, float *std, float *df, float proportion, float c0lim)
{
  int i, total=0, *dt, ntarget;
  float *n_t, *n_std, *n_df,  res;

  assert(dt=(int *)malloc(sizeof(int)*n));
 
  for(i=0;i<n;i++) {
    if(R_FINITE(t[i])) {
      dt[total]=i;
      total++;
     }
  }
  
  assert(n_t=(float *)malloc(total*sizeof(float)));
  assert(n_std=(float *)malloc(total*sizeof(float)));
  assert(n_df=(float *)malloc(total*sizeof(float)));
  for(i=0;i<total;i++) {
    n_t[i]=t[dt[i]];
    n_std[i]=std[dt[i]];
    n_df[i]=df[dt[i]];
  }
   
  ntarget=(int)ceil((double)proportion/2*total);
   Rprintf("%d ", ntarget);
   if(ntarget<1) res <- NA_REAL;
   else {
     float mc0=0, *ttop, *c1, *df1, *c0, qtarget;
     double *p0, *ptarget;
     int *index;

     assert(index=(int *)malloc(sizeof(int)*total));
     assert(ttop=(float *)malloc(sizeof(float)*ntarget));
     assert(c1=(float *)malloc(sizeof(float)*ntarget));
     assert(df1=(float *)malloc(sizeof(float)*ntarget));
     assert(p0=(double *)malloc(sizeof(double)*ntarget));
     assert(ptarget=(double *)malloc(sizeof(double)*ntarget));
     assert(c0=(float *)malloc(sizeof(float)*ntarget));
     memset(c0, 0, sizeof(float)*ntarget);

     for(i=0;i<total;i++) {
       n_t[i]=fabs(n_t[i]);
       index[i]=i;
     }

     order_index(n_t, index, total);
     qsort(n_t, total, sizeof(n_t[0]), distCompare);
        
     for(i=0;i<ntarget;i++) {
       int j=total-i-1;
       ttop[i]=n_t[j];
       
       c1[i]=n_std[index[j]]*n_std[index[j]];
       df1[i]=n_df[index[j]];
       p0[i]=pt(-(double)(ttop[i]), (double)(df1[i]), 1, 0);
       ptarget[i]=((i+0.5)/2/total-(1-proportion)*p0[i])/ proportion;
     
       if(ptarget[i]>p0[i]) {
	 qtarget=(float)qt((ptarget[i]), (double)(df1[i]), 1, 0);
	 c0[i]=c1[i]*(((ttop[i]/qtarget)*(ttop[i]/qtarget))-1);
       }
       if(c0[i]>c0lim) c0[i]=c0lim;
       mc0+=c0[i];
     }
     mc0/=(ntarget*1.0);
     res=mc0;
     free(ttop);
     free(c1);
     free(c0);
     free(df1);
     free(p0);
     free(ptarget);
     free(index);
     }

   free(n_t);
   free(n_std);
   free(n_df);
   free(dt);

   return(res);
}

void get_ebayes(double *d, int *pnrow, int *pncol, int *L, int *nL, float *T, float *B, float *proportion)
{
  GENE_DATA data;
  create_gene_data(d,pnrow,pncol,L,&data);
  /*print_gene_data(&data);*/
  compute_t_mod_stat(&data, L, T, nL);
  compute_t_mod_B(&data, L, B, proportion);
  free_gene_data(&data);
 
}

void get_B(double *d, int *pnrow, int *pncol, int *L, int *nL, float *B, float *proportion)
{
  GENE_DATA data;
  create_gene_data(d,pnrow,pncol,L,&data);
  compute_t_mod_B(&data, L, B, proportion);
  free_gene_data(&data);
 
}

void compute_t_mod_B(GENE_DATA *pdata, int *L, float *B, const void *extra)
{
  TMOD_DATA tmod;
  float s2_prior=0, df_prior=0, *s2_post, *df_total;
  float varpriorlim, varprior, *r, *kernel;
  float *T;
  float proportion=*(float *)extra;
  int i;
  int nrow=pdata->nrow;
  int nL=pdata->nL;
  
  assert(df_total=(float *)malloc(sizeof(int)*nrow));
  assert(s2_post=(float *)malloc(sizeof(float)*nrow));
  assert(r=(float *)malloc(sizeof(float)*nrow));
  assert(kernel=(float *)malloc(sizeof(float)*nrow));
  assert(T=(float *)malloc(sizeof(float)*nrow));
  
  create_tmod_data(&nrow, &tmod);
  if(nL==1) t1_mod_stat_func(pdata, L, &tmod);
  else t2_mod_stat_func(pdata, L, &tmod);
      
  fitFDist(tmod.sigma2, tmod.df_resid, nrow, &s2_prior, &df_prior);

  for(i=0;i<nrow;i++) {
    if((df_prior<FLT_MAX)&&(R_FINITE(tmod.df_resid[i]))) df_total[i]=tmod.df_resid[i]+df_prior;
    else if(df_prior==FLT_MAX)  df_total[i]=FLT_MAX;
    else df_total[i]=NA_REAL;
  
    if(df_prior==FLT_MAX) s2_post[i]=s2_prior;
    else if(tmod.df_resid[i]==0) s2_post[i]=(df_prior*s2_prior)/df_prior;
    else s2_post[i]=((float)tmod.df_resid[i]*tmod.sigma2[i]+df_prior*s2_prior)/((float)tmod.df_resid[i]+df_prior);
    if(!R_FINITE(tmod.stdev_unscale[i])) T[i]=NA_REAL;
    else T[i]=tmod.mean[i]/(tmod.stdev_unscale[i]*sqrt(s2_post[i]));
  }
  
  varpriorlim=10/s2_prior;
  varprior=tmixture(T, nrow, tmod.stdev_unscale, df_total, proportion, varpriorlim);
  if(!R_FINITE(varprior)) varprior=1/s2_prior;
  if(varprior<0.1/s2_prior) varprior=0.1/s2_prior;

  for (i=0;i<nrow;i++) {
    r[i]=(tmod.stdev_unscale[i]*tmod.stdev_unscale[i]+varprior)/(tmod.stdev_unscale[i]*tmod.stdev_unscale[i]);
    if(df_prior==FLT_MAX)   kernel[i]=T[i]*T[i]*(1-1/r[i])/2;
    else  kernel[i]=(1+df_total[i])/2*log((T[i]*T[i]+df_total[i])/(T[i]*T[i]/r[i]+df_total[i]));
    B[i]=log(proportion/(1-proportion))-log(r[i])/2+kernel[i];
    }
    
  free(s2_post);
  free(T);
  free(df_total);
  free(kernel);
  free(r);
  free_tmod_data(&tmod);
}
    
     
       
     
     
				     
    
    
    
    

  

  
