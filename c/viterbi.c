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
    int i, j, istate, jstate, isig;
    long t;
    double max_logprob, this_logprob;
    double altfreq_logprob, ratio_logprob, loggaussnorm;
    int max_index, index1, index2;
    double **eprob, **alpha, **beta, **gamma, mll;

    loggaussnorm = -log(2.0*3.14159*model_params->sigma_pi*model_params->sigma_ratio);

    /* Calculate logs of parameters for easier manipulation. */

    log_pi = alloc_double_vector(model_params->N);
    for (i = 0; i < model_params->N; i++) {
        log_pi[i] = log((model_params->pi)[i]);
    } 

    log_a = alloc_double_matrix(model_params->N, model_params->N);
    log_delta = alloc_double_matrix(model_params->T, 2*model_params->N);

    psi = alloc_int_matrix(model_params->T, 2*model_params->N);

    for (i = 0; i < model_params->N; i++) {
        for (j = 0; j < model_params->N; j++) {
            log_a[i][j] = log(model_params->a[i][j]);
        }
    }

    /* Initialization */
    for (i = 0; i < model_params->N; i++) {
        ratio_log_prob(model_params, observations[0], i, &ratio_logprob);
        log_delta[0][i] = log_pi[i] + ratio_logprob;
        log_delta[0][i + model_params->N] = log_delta[0][i];
        altfreq_log_prob(model_params, observations[0], i, -1, &altfreq_logprob);
        log_delta[0][i] += altfreq_logprob;
        altfreq_log_prob(model_params, observations[0], i, 1, &altfreq_logprob);
        log_delta[0][i + model_params->N] += altfreq_logprob;

        psi[0][i] = 0;
        psi[0][i + model_params->N] = 0;
        /* fprintf(stderr, "State %d has logprobs %lf/%lf\n", i, ratio_logprob, altfreq_logprob); */
    }

    /* Recursion */
    for (t = 1; t < model_params->T; t++) {
        for (i = 0; i < 2*model_params->N; i++) {
            istate = (i < model_params->N) ? i : i - model_params->N;
            isig = (i < model_params->N) ? -1 : 1;
            max_logprob = -99999999.0;
            for (j = 0; j < 2*model_params->N; j++) {
                jstate = (j < model_params->N) ? j : j - model_params->N;
                this_logprob = log_delta[t-1][j] + log_a[jstate][istate]; 
                if (this_logprob > max_logprob) {
                    max_logprob = this_logprob;
                    max_index = j;
                }
            }
            ratio_log_prob(model_params, observations[t], istate, &ratio_logprob);
            log_delta[t][i] = max_logprob + ratio_logprob;
            altfreq_log_prob(model_params, observations[t], istate, isig, &altfreq_logprob);
            log_delta[t][i] += altfreq_logprob;
            psi[t][i] = max_index;
        }
    }

    /* Termination */
    *log_prob = -99999999.0;
    state_indices[model_params->T - 1] = 1;
    for (i = 0; i < 2*model_params->N; i++) {
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
    free_int_matrix(psi, 2*model_params->N);

    fprintf(stderr, "Done with viterbi\n");

    /* Find probabilities by solving for gammas: */

    eprob = alloc_double_matrix(model_params->T, 2*model_params->N);
    alpha = alloc_double_matrix(model_params->T, 2*model_params->N);
    beta = alloc_double_matrix(model_params->T, 2*model_params->N);
    gamma = alloc_double_matrix(model_params->T, 2*model_params->N);

    calc_eprobs(model_params, observations, eprob);
    calc_alphas(model_params, observations, alpha, eprob, &mll);
    calc_betas(model_params, observations, beta, eprob);
    calc_gammas(model_params, alpha, beta, gamma, eprob);

    for (t = 0; t < model_params->T; t++) {
        index1 = state_indices[t];
        index2 = (index1 >= model_params->N) ? index1 - model_params->N : index1 + model_params->N;
        state_probabilities[t] = gamma[t][index1] + gamma[t][index2];
    }

    free_double_matrix(eprob, model_params->T);
    free_double_matrix(alpha, model_params->T);
    free_double_matrix(beta, model_params->T);
    free_double_matrix(gamma, model_params->T);
}
