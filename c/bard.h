/**************************************************************************
#
#  This software/database is "United States Government Work" under the terms of
#  the United States Copyright Act.  It was written as part of the authors'
#  official duties for the United States Government and thus cannot be
#  copyrighted.  This software/database is freely available to the public for
#  use without a copyright notice.  Restrictions cannot be placed on its present
#  or future use. 
#  
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and data, the National Human Genome Research
#  Institute (NHGRI) and the U.S. Government does not and cannot warrant the
#  performance or results that may be obtained by using this software or data.
#  NHGRI and the U.S.  Government disclaims all warranties as to performance,
#  merchantability or fitness for any particular purpose. 
#  
#  In any work or product derived from this material, proper attribution of the
#  authors as the source of the software or data should be made, using "NHGRI
#  Genome Technology Branch" as the citation. 
#
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>

typedef struct params {
    char *program;
    char *modelfile;
    char *obsfile;
    double min, max, maxratio, success;
    double minexparg;
    long nobins;
    int fixtrans;
    int derivatives;
    int verbose;
    char *fasta;
    char *bam;
    char *region;
    char *bedfile;
    int minqual;
    int mapqual;

} Params;

typedef struct statedata {
    int minor;
    int total;

} StateData;

typedef struct modelparams {
    int N;
    StateData *states;
    long T;
    double trans_prob;
    double **a;
    double *pi;
    double mu;
    double sigma;
    double rho_contam;

} ModelParams;

typedef struct observation {
    char *chr;
    long pos;
    long tumortotaldepth, normaltotaldepth, tumoraltdepth;

} Observation;

void read_model(char *filename, ModelParams **model_params);
void read_observations(char *filename, ModelParams **model_params, Observation **observations);
void filter_highcopy_observations(ModelParams **model_params, Observation **observations);
void run_viterbi(ModelParams *model_params, Observation *observations, int *state_indices, double *probs, double *log_prob);
void ratio_log_prob(ModelParams *model_params, Observation obs, int i, double *ans);
void write_states(int *states, double *probs, double log_prob, ModelParams *model_params, Observation *observations);
double run_baumwelch(ModelParams **model_params, Observation *observations);
void run_baumwelch_contam_optimization(ModelParams **model_params, Observation *observations);
void calc_bprobs(ModelParams *model_params, Observation *observations, double **bprob);
void calc_eprobs(ModelParams *model_params, Observation *observations, double **eprob);
void calc_alphas(ModelParams *model_params, Observation *observations, double **alpha, double **bprob, double **eprob, double *mll);
void calc_betas(ModelParams *model_params, Observation *observations, double **beta, double **bprob, double **eprob);
void calc_gammas(ModelParams *model_params, double **alpha, double **beta, double **gamma);
void calc_xis(ModelParams *model_params, double **alpha, double **beta, double ***xi, double **bprob, double **eprob);
void write_model(ModelParams *model_params);
void dbaum_dmuratio(ModelParams *model_params, Observation *observations, double **gamma, double *derivative);
void plot_mu_prob(ModelParams *model_params, Observation *observations);
void dbaum_dsigratio(ModelParams *model_params, Observation *observations, double **gamma, double *derivative);
void dbaum_dsigpi(ModelParams *model_params, Observation *observations, double **gamma, double *derivative);
void dbaum_dtransprob(ModelParams *model_params, Observation *observations, double ***xi, double *derivative);
void exp_log_prob(ModelParams *model_params, Observation *observations, double **gamma, double *elogp);
void copy_params(ModelParams *from_model_params, ModelParams **to_model_params);
void small_double_check(double *minexparg);
double ***alloc_double_tensor(int rows, int cols, int levels);
double **alloc_double_matrix(int rows, int cols);
double *alloc_double_vector(int entries);
int **alloc_int_matrix(int rows, int cols);
void free_double_tensor(double ***tensor, int rows, int cols);
void free_double_matrix(double **matrix, int rows);
void free_double_vector(double *vector);
void free_int_matrix(int **matrix, int rows);
int *alloc_int_vector(int entries);
void free_int_vector(int *vector);
double binprob(long successes, long trials, double probsuccess, double *pbin);
