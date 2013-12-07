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

#include "bard.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

extern Params *parameters;

void run_viterbi(ModelParams *model_params, Observation *observations, int *state_indices, double *state_probabilities, double *log_prob)
{
    double **log_delta;
    int **psi;
    double *log_pi;
    double **log_a;
    int i, j;
    long t;
    double max_logprob, this_logprob;
    int max_index, index;
    double **bprob, **eprob, **alpha, **beta, **gamma, mll;

    /* Find probabilities by solving for gammas: */
    fprintf(stderr, "Starting viterbi\n");

    bprob = alloc_double_matrix(model_params->T, model_params->N);
    eprob = alloc_double_matrix(model_params->T, model_params->N);
    alpha = alloc_double_matrix(model_params->T, model_params->N);
    beta = alloc_double_matrix(model_params->T, model_params->N);
    gamma = alloc_double_matrix(model_params->T, model_params->N);

    calc_bprobs(model_params, observations, bprob);
    calc_eprobs(model_params, observations, eprob);
    calc_alphas(model_params, observations, alpha, bprob, eprob, &mll);
    calc_betas(model_params, observations, beta, bprob, eprob);
    calc_gammas(model_params, alpha, beta, gamma);

    /* Calculate logs of parameters for easier manipulation. */

    log_pi = alloc_double_vector(model_params->N);
    for (i = 0; i < model_params->N; i++) {
        log_pi[i] = log((model_params->pi)[i]);
    } 

    log_a = alloc_double_matrix(model_params->N, model_params->N);
    for (i = 0; i < model_params->N; i++) {
        for (j = 0; j < model_params->N; j++) {
            log_a[i][j] = log(model_params->a[i][j]);
        }
    }

    log_delta = alloc_double_matrix(model_params->T, model_params->N);

    psi = alloc_int_matrix(model_params->T, model_params->N);

    /* Initialization */
    for (i = 0; i < model_params->N; i++) {
        log_delta[0][i] = log_pi[i] + log(bprob[0][i]) + log(eprob[0][i]);
        psi[0][i] = 0;
    }

    /* Recursion */
    for (t = 1; t < model_params->T; t++) {
        for (i = 0; i < model_params->N; i++) {
            max_logprob = -99999999.0;
            for (j = 0; j < model_params->N; j++) {
                this_logprob = log_delta[t-1][j] + log_a[j][i]; 
                if (this_logprob > max_logprob) {
                    max_logprob = this_logprob;
                    max_index = j;
                }
            }
            log_delta[t][i] = max_logprob + log(bprob[t][i]) + log(eprob[t][i]);
            psi[t][i] = max_index;
        }
    }

    /* Termination */
    *log_prob = -99999999.0;
    state_indices[model_params->T - 1] = 1;
    for (i = 0; i < model_params->N; i++) {
        if (log_delta[model_params->T - 1][i] > *log_prob) {
            *log_prob = log_delta[model_params->T - 1][i];
            state_indices[model_params->T - 1] = i;
        }
    }

    /* Traceback */
    for (t = model_params->T - 2; t >= 0; t--) {
        state_indices[t] = psi[t + 1][state_indices[t + 1]];
    }

    /* free memory */
    free_double_vector(log_pi);
    free_double_matrix(log_a, model_params->N);
    free_double_matrix(log_delta, model_params->T);
    free_int_matrix(psi, model_params->N);

    fprintf(stderr, "Done with viterbi\n");


    for (t = 0; t < model_params->T; t++) {
        index = state_indices[t];
        state_probabilities[t] = gamma[t][index];
        /* fprintf(stderr, "time %ld probability %lf for state %d\n", t, state_probabilities[t], index); */
    }

    free_double_matrix(eprob, model_params->T);
    free_double_matrix(alpha, model_params->T);
    free_double_matrix(beta, model_params->T);
    free_double_matrix(gamma, model_params->T);
}
