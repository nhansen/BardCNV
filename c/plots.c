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
    double **eprob; /* emission probabilities from each state with current parameters */
    double **alpha;
    double minratio, maxratio, mll;
   
    int i;
    long nopoints;

    minratio = (parameters->min != 0.0) ? parameters->min : 0.0;
    maxratio = (parameters->max != 0.0) ? parameters->max : 2.0;
    nopoints = parameters->nobins;

    eprob = (double **)malloc(sizeof(double *)*model_params->T);
    alpha = (double **)malloc(sizeof(double *)*model_params->T);

    for (i = 0; i < model_params->T; i++) {
        *(eprob + i) = (double *)malloc(sizeof(double)*2*model_params->N);
        *(alpha + i) = (double *)malloc(sizeof(double)*2*model_params->N);
    }

    fprintf(stderr, "mu from %lf to %lf, %d bins\n", minratio, maxratio, nopoints);
    for (i = 1; i < nopoints; i++) {
        model_params->mu_ratio = minratio + i*(maxratio - minratio)/(nopoints*1.0);
        fprintf(stderr, "Calculating probability for mu=%lf\n", model_params->mu_ratio);
        calc_eprobs(model_params, observations, eprob);
        calc_alphas(model_params, observations, alpha, eprob, &mll);
        fprintf(stdout, "%lf\t%lf\n", model_params->mu_ratio, mll);
    }
}

