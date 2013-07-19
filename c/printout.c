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

void write_states(int *states, double *probs, double log_prob, ModelParams *model_params, Observation *observations)
{
    int istate, qvalue;
    long t;
    double probability;
    StateData *thisState;
    Observation *obs;

    fprintf(stdout, "chr\tpos\tnormaltotaldepth\ttumortotaldepth\ttumoraltdepth\tminor\ttotal\tscore\n");
    for (t = 0; t < model_params->T; t++) {
        istate = *(states+t);
        probability = *(probs + t);
        /* fprintf(stdout, "prob %lf\n", probability); */
        qvalue = (probability == 1.0) ? 90 : -10 * ceil(log(1.0 - probability));
        if (qvalue > 90) {
            qvalue = 90;
        }
        thisState = model_params->states+istate;
        obs = &(observations[t]);
        fprintf(stdout, "%s\t%ld\t%ld\t%ld\t%ld\t%d\t%d\t%d\n", obs->chr, obs->pos, obs->normaltotaldepth, obs->tumortotaldepth, obs->tumoraltdepth, thisState->minor, thisState->total, qvalue);

    }
}

void write_model(ModelParams *model_params)
{
    int i, j;
    StateData *stateptr;
    double *pi;

    fprintf(stdout, "N= %d\n", model_params->N);
    stateptr = model_params->states;
    for (i=0; i<model_params->N; i++) {
        fprintf(stdout, "%d %d %g\n", stateptr->minor, stateptr->total, model_params->pi[i]);
        stateptr++;
    }

    /* fprintf(stdout, "Transitions\n");
    for (i=0; i<model_params->N; i++) {
        for (j=0; j<model_params->N; j++) {
            fprintf(stdout, "%lf", model_params->a[i][j]);
            if (j==model_params->N - 1) {
                fprintf(stdout, "\n");
            }
            else {
                fprintf(stdout, " ");
            }
        }
    } */
    fprintf(stdout, "Transprob= %g\n", model_params->trans_prob); 
    fprintf(stdout, "mu= %lf\n", model_params->mu); 
    fprintf(stdout, "sigma= %lf\n", model_params->sigma); 
    fprintf(stdout, "rho_contam= %lf\n", model_params->rho_contam); 
}
