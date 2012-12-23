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

/* "backward" procedure */

void calc_betas(ModelParams *model_params, Observation *observations, double **beta, double **eprob) {
    int i, istate, j, jstate, t, minor, total;
    double norm, sum, sig_ind_factor;
    double gaussnorm, exp_pi, exp_picontam, rho_contam;
    StateData *state;

    /* extra factor for normalization of Gaussian densities */
    gaussnorm = 1.0/(4.0*3.14159*model_params->sigma_pi*model_params->sigma_ratio);

    rho_contam = model_params->rho_contam;

    /* Initialization */
    for (i = 0; i < 2*model_params->N; i++) {
        beta[model_params->T - 1][i] = 0.5/model_params->N;
    }

    /* Induction */
    for (t = model_params->T - 2; t >= 0; t--) {
        norm = 0.0;
        for (i = 0; i < 2*model_params->N; i++) {
            istate = (i >= model_params->N) ? i - model_params->N : i;
            sum = 0.0;
            for (j = 0; j < model_params->N; j++) {
                sum += model_params->a[istate][j]/2.0 * (beta[t+1][j]*eprob[t+1][j]
                     + beta[t+1][j + model_params->N]*eprob[t+1][j + model_params->N]);
            } 
            beta[t][i] = sum;
            norm += sum;
        }
        for (i = 0; i < 2*model_params->N; i++) {
            beta[t][i] /= norm;
        }
    }
    /* fprintf(stderr, "Done with backward\n"); */
}

