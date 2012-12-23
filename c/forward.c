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

/* "forward" procedure */

void calc_alphas(ModelParams *model_params, Observation *observations, double **alpha, double **eprob, double *mll) {

    int i, j, t, minor, total;
    double norm, sum, gaussnorm, exp_pi, exp_picontam, rho_contam;
    StateData *state;

    /* extra factor for normalization of Gaussian densities */
    gaussnorm = 1.0/(4.0*3.14159*model_params->sigma_pi*model_params->sigma_ratio);

    rho_contam = model_params->rho_contam;
    
    /* Initialization */
    norm = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = &((model_params->states)[i]);
        minor = state->minor;
        total = state->total;

        exp_pi = (total==0) ? 0.5 : 1.0*minor/total;
        exp_picontam = rho_contam * (0.5 - exp_pi) + exp_pi;

        alpha[0][i] = model_params->pi[i] * eprob[0][i];
        alpha[0][i + model_params->N] = model_params->pi[i] * eprob[0][i + model_params->N];
        norm += alpha[0][i] + alpha[0][i + model_params->N];
        /* fprintf(stderr, "Initialization state %d eprob is %lf, %lf\n", i, eprob[0][i], eprob[0][i + model_params->N]); */
    }

    if (norm <= 0.0) {
        fprintf(stderr, "Underflow in alphas at t=0\n");
        exit(1);
    }

    for (i = 0; i < 2*model_params->N; i++) {
        alpha[0][i] /= norm;
    }

    *mll = log(norm);
    /* fprintf(stderr, "After initialization %lf\n", *mll); */

    /* Induction */
    for (t = 0; t < model_params->T - 1; t++) {
        norm = 0.0;
        for (j = 0; j < model_params->N; j++) {
            state = &((model_params->states)[j]);
            minor = state->minor;
            sum = 0.0;
            for (i = 0; i < model_params->N; i++) {
                sum += alpha[t][i] * model_params->a[i][j]/2.0; 
                sum += alpha[t][i + model_params->N] * model_params->a[i][j]/2.0; 
            }
            alpha[t + 1][j] = sum * eprob[t+1][j];
            alpha[t + 1][j + model_params->N] = sum * eprob[t+1][j + model_params->N];
            norm += alpha[t + 1][j] + alpha[t + 1][j + model_params->N];
        }
        if (norm <= 0.0) {
            fprintf(stderr, "Underflow in alphas at t=%d\n", t);
            exit(1);
        }
        sum = 0.0;
        for (j = 0; j < 2*model_params->N; j++) {
            alpha[t + 1][j] /= norm;
            sum += alpha[t + 1][j];
        }
        if (sum < 0.99 || sum > 1.01) {
            fprintf(stderr, "Ineffective normalization, resulting in sum %lf (norm %lf), at t=%d\n", sum, norm, t + 1);
        }
        if (norm <=0) {
            fprintf(stderr, "Zero norm for states at t=%d!\n", t);
            exit(1);
        }
        *mll += log(norm);
    }
    /* fprintf(stderr, "Done with forward\n"); */
}

