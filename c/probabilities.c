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

/* binomial probabilities */

double binprob(long successes, long trials, double probsuccess, double *pbin) {

    long i, high;
    double combprod, logcombprod, logpbin;

    if (trials - successes > successes) {
        high = successes;
    }
    else {
        high = trials - successes;
    }
    combprod = 1.0;
    logcombprod = 0.0;
    for (i = 0; i < high; i++) {
        combprod *= (double) (trials - i);
        combprod /= (double) (i + 1);
        combprod *= probsuccess * (1.0 - probsuccess);
        logcombprod += log(1.0*(trials - i)) - log(1.0*(i + 1)) + 
                     log(probsuccess) + log(1.0 - probsuccess);
    }

    /* *pbin = combprod * pow(probsuccess, successes-high) * pow(1.0-probsuccess, trials-successes-high); */

    logpbin = logcombprod + (successes - high) * log(probsuccess) + 
                     (trials - successes - high)*log(1.0-probsuccess);
    if (logpbin < parameters->minexparg) {
        *pbin = 0.0;
    }
    else {
        *pbin = exp(logpbin);
    }
    /* fprintf(stderr, "Probability of %ld out of %ld with psuccess %lf: %lf\n", successes, trials, probsuccess, *pbin); */
}

/* derivatives of Baum auxiliary function */

void dbaum_dmu(ModelParams *model_params, Observation *observations, double **gamma, double *derivative) {

    int i, total;
    long t;
    double sig, sum, addend, rho_contam, cneff;
    StateData *state;

    sig = model_params->sigma;
    rho_contam = model_params->rho_contam;
    sum = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        total = state->total;
        cneff = 2.0 * rho_contam + (1.0 - rho_contam) * total;

        for (t = 0; t < model_params->T; t++) {
            addend = gamma[t][i] * ((&(observations[t]))->tumortotaldepth -
                    cneff * model_params->mu * (&(observations[t]))->normaltotaldepth/2.0)/(sig*sig);
            sum += addend;
        }
    }
    *derivative = sum;
}

void dbaum_dsigma(ModelParams *model_params, Observation *observations, double **gamma, double *derivative) {

    int i, total;
    long t;
    double sig, sum, addend, thresh, rho_contam, cneff;
    StateData *state;

    thresh = 10.0;
    sig = model_params->sigma;
    rho_contam = model_params->rho_contam;
    sum = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        total = state->total;
        cneff = 2.0 * rho_contam + (1.0 - rho_contam) * total;
        for (t = 0; t < model_params->T; t++) {
            addend = gamma[t][i] * (-1.0/sig);
            sum += addend;
            addend = gamma[t][i] * 2.0 * pow((&(observations[t]))->tumortotaldepth - 
                         cneff * model_params->mu * (&(observations[t]))->normaltotaldepth/2.0, 2) / 
                            (pow(sig, 3) * cneff * (&(observations[t]))->normaltotaldepth);
            sum += addend;
        }
    }
    *derivative = sum;
}

void dbaum_dtransprob(ModelParams *model_params, Observation *observations, double ***xi, double *derivative) {

    int i, j;
    long t;
    double sum;

    sum = 0.0;
    for (t = 0; t < model_params->T - 1; t++) {
        for (i = 0; i < model_params->N; i++) {
            for (j = 0; j < model_params->N; j++) {
                if (i == j) {
                    sum += xi[t][i][j] * (1.0 - model_params->N) /
                          (1.0 - (model_params->N - 1.0) * model_params->trans_prob);
                }
                else {
                    sum += xi[t][i][j] / model_params->trans_prob;
                }
            }
        }
    }
    *derivative = sum;
}

/* WARNING: THIS ROUTINE IS NOT CORRECT SINCE CHANGEOVER FROM GAUSSIAN TO BINOMIAL ALT PROBABILITIES */
void exp_log_prob(ModelParams *model_params, Observation *observations, double **gamma, double *elogp) {

    int i, j, minor, total;
    long t;
    double sig, sum, addend, loggauss, rho_contam;
    StateData *state;

    sum = 0.0;
    rho_contam = model_params->rho_contam;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        minor = state->minor;
        total = state->total;

        loggauss = 0.5 * log(2.0*3.14159*model_params->sigma*model_params->sigma);

        for (t = 0; t < model_params->T; t++) {
            addend = gamma[t][i] * loggauss;
            sum -= addend;
        }
    }

    /* fprintf(stderr, "Done\n"); */
    *elogp = sum;
}
