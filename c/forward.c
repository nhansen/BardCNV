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

    int i, j, t;
    double norm, sum;

    /* Initialization */
    norm = 0.0;
    for (i = 0; i < model_params->N; i++) {

        alpha[0][i] = 0.5*model_params->pi[i] * eprob[0][i];
        alpha[0][i + model_params->N] = 0.5*model_params->pi[i] * eprob[0][i + model_params->N];
        norm += alpha[0][i] + alpha[0][i + model_params->N];
    }

    if (norm <= 0.0) {
        fprintf(stderr, "Underflow in alphas at t=0--assigning all probability to 0 state.\n");
        norm = 1.0;
        alpha[0][0] = 1.0;
    }

    for (i = 0; i < 2*model_params->N; i++) {
        alpha[0][i] /= norm;
    }

    *mll = log(norm);
    /* fprintf(stderr, "After initialization %lf\n", *mll); */

    /* Induction */
    for (t = 0; t < model_params->T - 1; t++) {
        norm = 0.0;
        for (j = 0; j < 2*model_params->N; j++) {
            sum = 0.0;
            for (i = 0; i < 2*model_params->N; i++) {
                sum += alpha[t][i] * model_params->a[i][j];
            }
            alpha[t + 1][j] = sum * eprob[t+1][j];
            norm += alpha[t + 1][j];
        }
        if (norm <= 0.0) {
            fprintf(stderr, "Underflow in alphas at t=%d--assigning all probability to 0 state.\n", t);
            norm = 1.0;
            alpha[t+1][0] = 1.0;
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

