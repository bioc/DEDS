#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <float.h>

#define EPSILON (12*FLT_EPSILON)

typedef struct tagGENE_DATA{
  float **d; /*the gene values matrix, mxn*/
  int nrow; /*nrow is the number of the genes*/
  int ncol; /*ncol is the number of the experiments*/
  int *L;   /*the status labelling of each experiment*/
  int nL;
} GENE_DATA;

typedef int (*FUNC_SAMPLE)(int *);
typedef void (*FUNC_COMPUTE_STAT)(GENE_DATA*, int*, float *, const void*);
typedef float(*FUNC_STAT)(const float *, const int *, const int, const void*);
typedef float (*FUNC_MAX)(float *, int);
typedef void (*FUNC_COMPUTE_P)(float *, float *, int *, int *, int *, int *, float *);

typedef struct tagTEST_DATA{
  int n_stat;
  int weighted_dist;
  float *extras;
  FUNC_COMPUTE_STAT *stat_array; 
  FUNC_COMPUTE_P func_compute_p;
  FUNC_MAX  func_max;
  FUNC_SAMPLE func_next_sample;
} TEST_DATA;

typedef struct tagDEDS_RES{
  int nsig;
  int nT;
  int nrow;
  int *R;
  float *E;
  float *D;
  double *FDR;
  float **T;
  float *wval;
} DEDS_RES; 

typedef struct tagTMOD_DATA{
  int nrow;
  float *mean;
  float *sigma2;
  int *df_resid;
  float *stdev_unscale;
} TMOD_DATA; 

/*********************************************************************/
/*               utilities                                           */
/*********************************************************************/
float max_high(float *X, int n);
float max_low(float *X, int n);
float max_abs(float *X, int n);
void compute_euclid(float **X, int nrow, int ncol, float *E, float *wval, float *dist);
void order_index(float *V,int *R,int n);
void order_dist(float *V,int n);
int indexCompare(const void *v1, const void *v2) ;
int distCompare(const void *v1, const void *v2) ;
float select(unsigned long k, unsigned long n, float arr[]);
float median(float *X, int n);
void quantile(float *X, int nX, float *q, int nq, float *ret);
float dpowern(float n, int k);
float mad(float *X, int n);
float trigammaInverse(float x);
void fitFDist(float *sigma2, int *df1, int n, float *scale, float *df2);
void inv_chol(double *x, int *nr, int *nc, double *v, int *info);
float tmixture(float *t, int n, float *std, float *df, float proportion, float c0lim);
/*********************************************************************/
/*               format input data                                   */
/*********************************************************************/
void malloc_gene_data(GENE_DATA *pdata);
void free_gene_data(GENE_DATA *pdata);
void create_gene_data(double*d,  int*pnrow, int *pncol, int *L, GENE_DATA *pdata);
void print_gene_data(GENE_DATA *pdata);
void print_b(int b, int B, char *prompt);
void create_deds_res(int *pnrow, int *pnsig, int *pnT, DEDS_RES *pdr);
void free_deds_res(DEDS_RES *pdr);
void extract_deds_res(DEDS_RES *pdr, double *E, int *R, double *FDR, double *T);
void create_tmod_data(int *pnrow, TMOD_DATA *ptmod);
void free_tmod_data(TMOD_DATA *ptmod);
/*********************************************************************/
/*               link R                                              */
/*********************************************************************/
FUNC_COMPUTE_STAT type2stat(char *ptest, int *nL);
int type2test(char **options, TEST_DATA *td, int *nT, int *nL, float *extras);
void get_stat(double *d, int *pnrow, int *pncol, int *L, float *T,char **options, float *extras, int *nL);
void get_unadjp(double *d, int *pnrow, int *pncol, int *L, float *T, float *P, char **options, float *extras, int *nL, int *B);
void get_deds_order(double *d, int *pnrow, int *pncol, int *L, char **options, float *extras, int *nL, int *nT, int *B, double *E, int *R, double *T);
void func_get_order(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B);
void func_get_FDR(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B);
void func_deds_quick(GENE_DATA *pdata, TEST_DATA *ptd, DEDS_RES *pdr, int *B) ;
void get_deds_FDR(double *d, int *pnrow, int *pncol, int *L, char **options, float *extras, int *quick, 
		   int *nL, int *nT, int *B, int *pnsig, double *E, int *R, double *FDR, double *T);
void calc_FDR(float *bD, float *D, int *R, int *pnrow, int *pncol, int *nsig, float *F);
void calc_adjP(float *bD, float *D, int *R, int *pnrow, int *pncol, int *nsig, float *F);
void get_ebayes(double *d, int *pnrow, int *pncol, int *L, int *nL, float *T, float *B, float *proportion);
void get_B(double *d, int *pnrow, int *pncol, int *L, int *nL, float *B, float *proportion);
void get_t_mod_stat(double *d, int *pnrow, int *pncol, int *L, float *T, int *nL);
void get_f_mod_stat(double *d, int *pnrow, int *pncol, int *L, float *T, int *nL);

/*********************************************************************/
/*                             tests                                 */
/*********************************************************************/
float t2_stat(const float *Y, const int *L, const int n, const void *extra);
float t1_stat(const float *Y, const int *L, const int n, const void *extra);
float fc2_stat(const float *Y, const int *L, const int n, const void *extra);
float fc1_stat(const float *Y, const int *L, const int n, const void *extra);
float fcm_stat(const float *Y, const int *L, const int n,  const void *extra);
float f_stat(const float *Y, const int *L, const int n, const void *extra);
void compute_t2_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_t1_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_fc2_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_fc1_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_fcm_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_sam2_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_sam1_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_samm_stat_q(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_f_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_t_mod_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_f_mod_stat(GENE_DATA *pdata, int *L, float *T, const void *extra);
void compute_t_mod_B(GENE_DATA *pdata, int *L, float *B, const void *extra);
void t1_mod_stat_func(GENE_DATA *pdata, int *L, TMOD_DATA *ptmod);
void t2_mod_stat_func(GENE_DATA *pdata, int *L, TMOD_DATA *ptmod);

/*********************************************************************/
/*                             Sampling                              */
/*********************************************************************/
int bincoeff(int n, int k);
void set_seed(long int seed);
float get_rand();
void sample(int *V, int n, int m);
int first_sample(int *L);
int next_sample(int *L);
int next_sample_1(int *L);
void sample2label(int n, int k, int *nk,int *permun, int *L);
void creat_sampling(int n, int *L, int B);
void delete_sampling();

