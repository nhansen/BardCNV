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

void plot_mu_prob(ModelParams *model_params, Observation *observations)
{
    double **bprob; /* emission probabilities from each state with current parameters */
    double **eprob; /* emission probabilities from each state with current parameters */
    double **alpha;
    double minmu, maxmu, mll;
   
    int i;
    long nopoints;

    minmu = (parameters->min != 0.0) ? parameters->min : 0.0;
    maxmu = (parameters->max != 0.0) ? parameters->max : 2.0;
    nopoints = parameters->nobins;

    bprob = alloc_double_matrix(model_params->T, model_params->N);
    eprob = alloc_double_matrix(model_params->T, model_params->N);
    alpha = alloc_double_matrix(model_params->T, model_params->N);

    fprintf(stderr, "mu from %lf to %lf, %d bins\n", minmu, maxmu, nopoints);
    for (i = 1; i < nopoints; i++) {
        model_params->mu = minmu + i*(maxmu - minmu)/(nopoints*1.0);
        fprintf(stderr, "Calculating probability for mu=%lf\n", model_params->mu);
        calc_bprobs(model_params, observations, bprob);
        calc_eprobs(model_params, observations, eprob);
        calc_alphas(model_params, observations, alpha, bprob, eprob, &mll);
        fprintf(stdout, "%lf\t%lf\n", model_params->mu, mll);
    }

    free_double_matrix(eprob, model_params->T);
    free_double_matrix(alpha, model_params->T);
}
