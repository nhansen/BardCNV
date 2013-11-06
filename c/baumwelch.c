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

#define TOL 0.000001

extern Params *parameters;

double run_baumwelch(ModelParams **model_params, Observation *observations)
{
    double **bprob; /* emission probabilities of B-allele freq from each state (param indep) */
    double **eprob; /* emission probabilities for tumor read depth from each state with current parameters */
    double **alpha; /* "forward" variables (see Rabiner) */
    double **beta; /* "backward" variables (see Rabiner) */
    double **gamma; /* single-site state probabilities */
    double ***xi; /* two-site state probabilities */
    double initial_mll, final_mll, mll_diff, elogp, elogpnew;
    double norm, sum, mu_num, mu_den, random_number, sigma_num, small_sigma_num, sigma_add;
    double addend, cneff;
    double derivative, last_mu, last_sigma, rho_contam, new_diag;
   
    StateData *thisstate; 
    int i, j, total, bt_iterations;
    long t;

    rho_contam = (*model_params)->rho_contam;

    bprob = alloc_double_matrix((*model_params)->T, (*model_params)->N);
    eprob = alloc_double_matrix((*model_params)->T, (*model_params)->N);
    alpha = alloc_double_matrix((*model_params)->T, (*model_params)->N);
    beta = alloc_double_matrix((*model_params)->T, (*model_params)->N);
    gamma = alloc_double_matrix((*model_params)->T, (*model_params)->N);

    if (parameters->fixtrans == 0) {
        xi = alloc_double_tensor((*model_params)->T, (*model_params)->N,
                             (*model_params)->N);
    }

    calc_bprobs(*model_params, observations, bprob); /* binary emission probs--only need to calculate once */
    calc_eprobs(*model_params, observations, eprob); /* will be recalculated with each iteration of BW */
    calc_alphas(*model_params, observations, alpha, bprob, eprob, &final_mll);
    calc_betas(*model_params, observations, beta, bprob, eprob);
    initial_mll = (1.0 - 2.0*TOL) * final_mll; /* to assure at least one iteration */
    calc_gammas(*model_params, alpha, beta, gamma);
    if (parameters->fixtrans == 0) {
        calc_xis(*model_params, alpha, beta, xi, bprob, eprob);
    }

    mll_diff = final_mll - initial_mll;
    bt_iterations = 0; 

    /* pre-calculate numerator of most likely mu (of tumor read depth), */
    /* since it is independent of state probabilities: */

    mu_num = 0.0;
    for (t = 0; t < (*model_params)->T; t++) {
        mu_num += (double)(&(observations[t]))->tumortotaldepth;
    }

    fprintf(stderr, "Initial/final mll %lf %lf\n", initial_mll, final_mll);

    while ((mll_diff)/final_mll > TOL || (-1.0*mll_diff)/final_mll > TOL) {

        initial_mll = final_mll;
      
        /* store last */
        last_mu = (*model_params)->mu; 
        last_sigma = (*model_params)->sigma; 

        /* expected value of log(P) */
        /* exp_log_prob(*model_params, observations, gamma, &elogp); */

        /* re-estimate */
        /* initial probabilities */
        for (i = 0; i < (*model_params)->N; i++) {
            (*model_params)->pi[i] = gamma[0][i];
            
        }

        /* transition probabilities and sums for emission parameters */

        mu_den = 0.0;

        /* first take care of transition probabilities */
        norm = 0.0;
        sum = 0.0;
        for (i = 0; i < (*model_params)->N; i++) {
            if (parameters->fixtrans == 0) {
                for (j = 0; j < (*model_params)->N; j++) {
                    for (t = 0; t < (*model_params)->T - 1; t++) {
                        addend = xi[t][i][j];
                        if (j != i) {
                            sum += addend;
                        }
                        norm += addend;
                    }
                }
            }

            /* state info for estimating ratio and binomial parameters */
            thisstate = &(((*model_params)->states)[i]);
            total = thisstate->total;
            cneff = 2.0 * rho_contam + (1.0 - rho_contam) * total;
            
            for (t = 0; t < (*model_params)->T; t++) {
                mu_den += gamma[t][i] * (cneff * 0.5)
                                       * (&(observations[t]))->normaltotaldepth;
            }
        }
        if (parameters->fixtrans == 0) {
            (*model_params)->trans_prob = sum / (norm * (((*model_params)->N) - 1.0));
            fprintf(stderr, "New trans_prob value of %lf\n", (*model_params)->trans_prob);
            for (i = 0; i < (*model_params)->N; i++) {
                for (j = 0; j < (*model_params)->N; j++) {
                    if (i != j) {
                        ((*model_params)->a)[i][j] = (*model_params)->trans_prob;
                    }
                    else {
                        new_diag = 1.0 - (*model_params)->trans_prob * 
                                          ((*model_params)->N - 1.0);
                        ((*model_params)->a)[i][j] = new_diag;
                    }
                }
            }
        }
        (*model_params)->mu = mu_num/mu_den;
        fprintf(stderr, "New mu value of %lf (num %lf, den %lf)\n", (*model_params)->mu, mu_num, mu_den);

        sigma_num = 0.0;
        small_sigma_num = 0.0;
        for (i = 0; i < (*model_params)->N; i++) {

            thisstate = &(((*model_params)->states)[i]);
            total = thisstate->total;

            cneff = 2.0 * rho_contam + (1.0 - rho_contam) * total;
           
            for (t = 0; t < (*model_params)->T; t++) {
                addend = pow( ( (&(observations[t]))->tumortotaldepth -
                    cneff*(&(observations[t]))->normaltotaldepth*(*model_params)->mu/2.0 ), 2) *
                              gamma[t][i] / ( cneff * (&(observations[t]))->normaltotaldepth );
                if (sigma_num == 0.0 || addend/sigma_num > 0.000000001) {
                    sigma_num += addend;
                }
                else { /* aggregate small terms */
                    small_sigma_num += addend;
                }
            }
        }

        (*model_params)->sigma = sqrt(2.0 * (sigma_num + small_sigma_num) / (*model_params)->T);
        fprintf(stderr, "New sigma value of %lf\n", (*model_params)->sigma);

        if (parameters->derivatives != 0) {
            /* This section is an optional derivative check to assure precision is still good */
            dbaum_dmu(*model_params, observations, gamma, &derivative);
            fprintf(stderr, "Mu derivative is %lf\n", derivative);
            dbaum_dsigma(*model_params, observations, gamma, &derivative);
            fprintf(stderr, "Sig derivative is %lf\n", derivative);
            if (parameters->fixtrans == 0) {
                dbaum_dtransprob(*model_params, observations, xi, &derivative);
                fprintf(stderr, "Transprob derivative is %lf\n", derivative);
            }
    
            /* exp_log_prob(*model_params, observations, gamma, &elogpnew); */
            /* if (elogpnew < elogp) {
                fprintf(stderr, "LOWER EXPECTED PROB!!!!\n");
            }
            fprintf(stderr, "New/Old expected logp values with old states: %lf, %lf\n", elogp, elogpnew); */
        }

        calc_eprobs(*model_params, observations, eprob);
        calc_alphas(*model_params, observations, alpha, bprob, eprob, &final_mll);
        calc_betas(*model_params, observations, beta, bprob, eprob);
        calc_gammas(*model_params, alpha, beta, gamma);
        if (parameters->fixtrans == 0) {
            calc_xis(*model_params, alpha, beta, xi, bprob, eprob);
        }

        mll_diff = final_mll - initial_mll;
        if (mll_diff >= 0) {
            fprintf(stderr, "MLL difference: %lf (%lf, %lf)\n", mll_diff, initial_mll, final_mll);
        }
        else {
            while (mll_diff < 0) {
                fprintf(stderr, "ERROR!!!BACKTRACKING: MLL difference: %lf (%lf, %lf)\n", mll_diff, initial_mll, final_mll);
                random_number = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
                (*model_params)->mu = last_mu + random_number;
                random_number = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
                (*model_params)->sigma = last_sigma + random_number;
                calc_eprobs(*model_params, observations, eprob);
                calc_alphas(*model_params, observations, alpha, bprob, eprob, &final_mll);
                calc_betas(*model_params, observations, beta, bprob, eprob);
                calc_gammas(*model_params, alpha, beta, gamma);
                if (parameters->fixtrans == 0) {
                    calc_xis(*model_params, alpha, beta, xi, bprob, eprob);
                }
                mll_diff = final_mll * TOL + 1;
                bt_iterations++;
                if (bt_iterations > 3) {
                    break;
                }
            }
        }
        if (bt_iterations > 3) {
            break;
        }
    }

    /* free memory */
    free_double_matrix(bprob, (*model_params)->T);
    free_double_matrix(eprob, (*model_params)->T);
    free_double_matrix(alpha, (*model_params)->T);
    free_double_matrix(beta, (*model_params)->T);
    free_double_matrix(gamma, (*model_params)->T);
    if (parameters->fixtrans == 0) {
        free_double_tensor(xi, (*model_params)->T, (*model_params)->N);
    }

    return final_mll;
}

/* pre-calculate the probability of the observed number of alternate allele reads, given the total tumor read depth from a particular state at a particular time */

void calc_bprobs(ModelParams *model_params, Observation *observations, double **bprob) {
    int i, total, minor;
    long t;
    long tumordepth, altdepth;
    double rho_contam, abinaryprob, bbinaryprob;
    double *propeff; /* expected proportion of B-allele, given contamination, for a state */
    StateData *state;

    rho_contam = model_params->rho_contam;
   
    propeff = alloc_double_vector(model_params->N);

    /* pre-calculate expected ratios for each state */ 
    for (i = 0; i < model_params->N; i++) {
        state = &((model_params->states)[i]);
        total = state->total;
        minor = state->minor;
        propeff[i] = (rho_contam + (1.0-rho_contam)*minor)/(2.0*rho_contam + (1.0-rho_contam)*total);

        for (t=0; t < model_params->T; t++) {
            tumordepth = (&(observations[t]))->tumortotaldepth;
            altdepth = (&(observations[t]))->tumoraltdepth;
            
            binprob(altdepth, tumordepth, propeff[i], &bbinaryprob);
            binprob(altdepth, tumordepth, 1.0 - propeff[i], &abinaryprob);

            /* fprintf(stderr, "Calculated probs %lf %lf for depth %d out of %d with prob %lf\n", bbinaryprob, abinaryprob, altdepth, tumordepth, propeff[i]); */
            bprob[t][i] = 0.5*(abinaryprob + bbinaryprob);
            if (isnan(bprob[t][i])) {
                fprintf(stderr, "bprob is not a number at t=%ld, state %d!\n", t, i);
                exit(1);
            }
            if (bprob[t][i] <= 0.00001) {
                if ( parameters->verbose != 0 ) {
                    fprintf(stderr, "Setting bprob to 0.00001 at %ld for state %d, tumoralt %ld, tumortotal %ld (prop=%lf)\n", t, i, altdepth, tumordepth, propeff[i]);
                }
                bprob[t][i] = 0.00001;

            }
            if (bprob[t][i] > 1.0) {
                fprintf(stderr, "bprob is greater than 1.0 at t=%ld, state %d (tumor %ld, normal %ld, prop %lf) probability; %lf!\n", t, i, tumordepth, altdepth, propeff[i], bprob[t][i]);
                exit(1);
            }
        }

        /* fprintf(stderr, "Finished setting binary probabilities for state %d.\n", i); */
    }

    /* check for probable states */

    /* for (t=0; t < model_params->T; t++) {
        bbinaryprob = 0;
        for (i = 0; i < model_params->N; i++) {
            if (bprob[t][i] > bbinaryprob) {
                bbinaryprob = bprob[t][i];
            }
        }

        fprintf(stderr, "Highest binary probability %lf at t=%d\n", bbinaryprob, t);
    } */
}

/* pre-calculate the probability of the tumor read depth observation from a particular state at a particular time */

void calc_eprobs(ModelParams *model_params, Observation *observations, double **eprob) {
    int i, total;
    long t;
    long tumordepth, normaldepth;
    double mu, sigma, rho_contam;
    double gaussnorm, logprob;
    double *statenorm, *cneff;
    StateData *state;

    mu = model_params->mu;
    sigma = model_params->sigma;
    rho_contam = model_params->rho_contam;

    /* extra factor for normalization of Gaussian densities */    
    gaussnorm = 1.0/sqrt(3.14159)/sigma; /* will include cneff and normaltotaldepth later */

    /* pre-calculate normalization factor and effective copy number for each different state: */
    statenorm = alloc_double_vector(model_params->N);
    cneff = alloc_double_vector(model_params->N);

    for (i = 0; i < model_params->N; i++) {
        state = &((model_params->states)[i]);
        total = state->total;
        cneff[i] = 2.0*rho_contam + (1.0-rho_contam)*total;

        statenorm[i] = gaussnorm/sqrt(cneff[i]);
    }

    for (t = 0; t < model_params->T; t++) {
        tumordepth = (&(observations[t]))->tumortotaldepth;
        normaldepth = (&(observations[t]))->normaltotaldepth;
        for (i = 0; i < model_params->N; i++) {
            logprob = pow((1.0*tumordepth - mu*cneff[i]*normaldepth/2.0), 2);
            logprob = -1.0*logprob/(normaldepth*cneff[i]*pow(sigma, 2));
            if (logprob < parameters->minexparg) {
                if ( parameters->verbose != 0 ) {
                    fprintf(stderr, "Setting eprob to 0.0 at %ld for state %d, normal %ld, tumor %ld (mu=%lf, sigma=%lf, cneff=%lf)\n", t, i, normaldepth, tumordepth, mu, sigma, cneff[i]);
                }
                eprob[t][i] = 0.0;
            }
            else {
                eprob[t][i] = exp(logprob) * statenorm[i] / sqrt(normaldepth);
            }
        }
    }
    free_double_vector(statenorm);
    free_double_vector(cneff);
    /* fprintf(stderr, "Done with eprobs\n"); */
}

void calc_gammas(ModelParams *model_params, double **alpha, double **beta, double **gamma) {
    int i, j;
    long t;
    double norm;

    for (t = 0; t < model_params->T; t++) {
        norm = 0.0;
        for (j = 0; j < model_params->N; j++) {
            gamma[t][j] = alpha[t][j] * beta[t][j];
            norm += gamma[t][j];
        }
        for (j = 0; j < model_params->N; j++) {
            gamma[t][j] /= norm;
        }
    }
}

void calc_xis(ModelParams *model_params, double **alpha, double **beta, double ***xi, double **bprob, double **eprob) {
    int i, j;
    long t;
    double norm;

    for (t = 0; t < model_params->T - 1; t++) {
        norm = 0.0;
        for (i = 0; i < model_params->N; i++) {
            for (j = 0; j < model_params->N; j++) {
                xi[t][i][j] = alpha[t][i] * beta[t+1][j] * model_params->a[i][j] * bprob[t+1][j] * eprob[t+1][j];
                norm += xi[t][i][j];
            }
        }
        for (i = 0; i < model_params->N; i++) {
            for (j = 0; j < model_params->N; j++) {
                xi[t][i][j] /= norm;
                /* fprintf(stdout, "XI: %d %d %lf\n", i, j, xi[t][i][j]); */
            }
        }
    }
    /* fprintf(stderr, "Done with xi\n"); */
}

void run_baumwelch_contam_optimization(ModelParams **model_params, Observation *observations)
{
    ModelParams *bestmodel;
    double min_contam, max_contam, contam, mll, best_mll;

    min_contam = parameters->min;
    max_contam = parameters->max;
    best_mll = -9999999999;
    for (contam = min_contam; contam <= max_contam; contam += 0.01) {
        (*model_params)->rho_contam = contam;
        fprintf(stderr, "Optimizing for contam %lf.\n", (*model_params)->rho_contam);
        mll = run_baumwelch(model_params, observations);
        fprintf(stderr, "MLL at contam %lf is %lf\n", contam, mll);
        if (mll > best_mll) {
            best_mll = mll;
            copy_params(*model_params, &bestmodel);
        }
    }

    copy_params(bestmodel, model_params);
    fprintf(stderr, "Optimal MLL at contam %lf is %lf.  Outputting model\n", (*model_params)->rho_contam, best_mll);
}

