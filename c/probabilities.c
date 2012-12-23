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

void ratio_log_prob(ModelParams *model_params, Observation obs, int i, double *ans)
{
    StateData *state;
    int total;
    double depthratio;
    double mu_ratio, sigma_ratio, rho_contam;
    double log_prob, diff;

    depthratio = (&obs)->depthratio;

    state = &((model_params->states)[i]);
    total = state->total;

    mu_ratio = model_params->mu_ratio;
    sigma_ratio = model_params->sigma_ratio;
    rho_contam = model_params->rho_contam;

    log_prob = 0.0;
    diff = depthratio - ( (1.0 - rho_contam)*total/2.0 + rho_contam ) * mu_ratio;
    log_prob -= 0.5*(diff*diff)/((rho_contam + (1.0 - rho_contam)*total/2.0) * 
                            sigma_ratio*sigma_ratio);

    /*fprintf(stderr, "State %d, ratio %lf logprob %lf\n", i, depthratio, log_prob); */
    *ans = log_prob;
}

void altfreq_log_prob(ModelParams *model_params, Observation obs, int i, int sig_pm, double *ans)
{
    StateData *state;
    int minor, total;
    double pi;
    double log_prob, diff;
    double sigma_pi, rho_contam;
    double exp_pi, exp_picontam;

    state = &((model_params->states)[i]);
    minor = state->minor;
    total = state->total;

    sigma_pi = model_params->sigma_pi;
    rho_contam = model_params->rho_contam;

    pi = (&obs)->pi;

    if (total > 0) {
        exp_pi = (sig_pm == -1) ? 1.0*minor/total : 1.0 - 1.0*minor/total;
        exp_picontam = rho_contam * (0.5 - exp_pi) + exp_pi;
    }
    else {
        exp_picontam = 0.5;
    }

    diff = pi - exp_picontam;
    log_prob = -1.0*(diff*diff)/(8.0*exp_picontam*(1.0-exp_picontam)*sigma_pi*sigma_pi);

    /* fprintf(stderr, "State %d (minor %d, total %d), sigpm %d, ratio %lf, diff %lf, sigma_pi %lf, logprob %lf\n", i, minor, total, sig_pm, pi, diff, sigma_pi, log_prob); */
    *ans = log_prob;
}

/* derivatives of Baum auxiliary function */

void dbaum_dmuratio(ModelParams *model_params, Observation *observations, double **gamma, double *derivative) {

    int i, j, t, total;
    double sig, sum, addend, rho_contam, eff_cn;
    StateData *state;

    sig = model_params->sigma_ratio;
    rho_contam = model_params->rho_contam;
    sum = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        total = state->total;
        eff_cn = rho_contam + (1.0 - rho_contam) * total/2.0;

        for (t = 0; t < model_params->T; t++) {
            addend = (gamma[t][i] + gamma[t][i + model_params->N]) *
                    ((&(observations[t]))->depthratio -
                    eff_cn * model_params->mu_ratio)/(sig*sig);
            /* if (addend > 0.01) {
                fprintf(stderr, "Large term t=%d, i=%d, %lf\n", t, i, addend);
            } */
            sum += addend;
        }
    }
    *derivative = sum;
}

void dbaum_dsigratio(ModelParams *model_params, Observation *observations, double **gamma, double *derivative) {

    int i, j, t, total;
    double sig, sum, addend, thresh, rho_contam, eff_cn;
    StateData *state;

    thresh = 10.0;
    sig = model_params->sigma_ratio;
    rho_contam = model_params->rho_contam;
    sum = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        total = state->total;
        eff_cn = rho_contam + (1.0 - rho_contam) * total/2.0;
        for (t = 0; t < model_params->T; t++) {
            addend = (gamma[t][i] + gamma[t][i + model_params->N]) * (-1.0/sig);
            sum += addend;
            addend = (gamma[t][i] + gamma[t][i + model_params->N]) *
                      pow((&(observations[t]))->depthratio - 
                         eff_cn * model_params->mu_ratio, 2) / (pow(sig, 3) * eff_cn);
            /* if (addend > 0.01) {
                fprintf(stderr, "Large term t=%d, i=%d, %lf\n", t, i, addend);
            } */
            sum += addend;
            /* if (sum > thresh || sum < -1.0*thresh) {
                fprintf(stderr, "Sum larger than %lf at t=%d, i=%d, %lf\n", thresh, t, i, sum);
                thresh *= 10.0;
            } */
        }
    }
    *derivative = sum;
}

void dbaum_dsigpi(ModelParams *model_params, Observation *observations, double **gamma, double *derivative) {

    int i, j, t, minor, total;
    double sig, sum, addend, rho_contam, exp_pi, exp_picontam;
    StateData *state;

    sig = model_params->sigma_pi;
    rho_contam = model_params->rho_contam;
    sum = 0.0;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        minor = state->minor;
        total = state->total;

        exp_pi = (total == 0) ? 0.5 : 1.0*minor/total;
        exp_picontam = rho_contam * (0.5 - exp_pi) + exp_pi;

        for (t = 0; t < model_params->T; t++) {
            addend = (gamma[t][i] + gamma[t][i + model_params->N]) * (-1.0/sig);
            sum += addend;
            addend = gamma[t][i] *
                      pow((&(observations[t]))->pi - exp_picontam, 2) 
                      * 0.25 / (pow(sig, 3) * exp_picontam * (1.0 - exp_picontam));
            sum += addend;
            addend = gamma[t][i + model_params->N] *
                      pow((&(observations[t]))->pi - (1.0 - exp_picontam), 2) 
                      * 0.25 / (pow(sig, 3) * exp_picontam * (1.0 - exp_picontam));
            sum += addend;
        }
    }
    *derivative = sum;
}

void exp_log_prob(ModelParams *model_params, Observation *observations, double **gamma, double *elogp) {

    int i, j, t, minor, total;
    double sig, sum, addend, logrp, logpip, logpin, loggauss, rho_contam, exp_pi, exp_picontam;
    StateData *state;

    sum = 0.0;
    rho_contam = model_params->rho_contam;
    for (i = 0; i < model_params->N; i++) {
        state = model_params->states + i;
        minor = state->minor;
        total = state->total;
        exp_pi = (total == 0) ? 0.5 : 1.0*minor/total;
        exp_picontam = rho_contam * (0.5 - exp_pi) + exp_pi;

        loggauss = log(4.0*3.14159*model_params->sigma_pi *
                                      model_params->sigma_ratio);
        loggauss += 0.5 * log(2.0*rho_contam + (1.0 - rho_contam) * total);
        loggauss += 0.5 * log(exp_picontam * (1.0 - exp_picontam));

        for (t = 0; t < model_params->T; t++) {
            ratio_log_prob(model_params, observations[t], i, &logrp);
            altfreq_log_prob(model_params, observations[t], i, -1, &logpin);
            altfreq_log_prob(model_params, observations[t], i, 1, &logpip);
            addend = gamma[t][i] * (logrp + logpin);
            sum += addend;
            addend = gamma[t][i + model_params->N] * (logrp + logpip);
            sum += addend;
            addend = (gamma[t][i] + gamma[t][i + model_params->N]) * loggauss;
            sum -= addend;
        }
    }

    *elogp = sum;
}
