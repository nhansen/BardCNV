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
    double **eprob; /* emission probabilities from each state with current parameters */
    double **alpha; /* "forward" variables (see Rabiner) */
    double **beta; /* "backward" variables (see Rabiner) */
    double **gamma; /* single-site state probabilities */
    double ***xi; /* two-site state probabilities */
    double initial_mll, final_mll, mll_diff, elogp, elogpnew;
    double norm, sum, mu_ratio_num, mu_ratio_den, random_number;
    double explow, exphigh, difflow, diffhigh;
    double total_pialt_num, total_pialt_den, addend, eff_cn;
    double sig_ratio_num, sig_ratio_den, small_sig_ratio_num, small_sig_ratio_den;
    double derivative, last_mu_ratio, last_sig_ratio, last_sig_pi, rho_contam, new_diag;
   
    StateData *thisstate; 
    int i, j, istate, jstate, minor, total, bt_iterations;
    long t;

    rho_contam = (*model_params)->rho_contam;

    eprob = alloc_double_matrix((*model_params)->T, 2*(*model_params)->N);
    alpha = alloc_double_matrix((*model_params)->T, 2*(*model_params)->N);
    beta = alloc_double_matrix((*model_params)->T, 2*(*model_params)->N);
    gamma = alloc_double_matrix((*model_params)->T, 2*(*model_params)->N);

    if (parameters->fixtrans == 0) {
        xi = alloc_double_tensor((*model_params)->T, 2*(*model_params)->N,
                             2*(*model_params)->N);
    }

    calc_eprobs(*model_params, observations, eprob);
    calc_alphas(*model_params, observations, alpha, eprob, &final_mll);
    initial_mll = (1.0 - 2.0*TOL) * final_mll; /* to assure at least one iteration */
    calc_betas(*model_params, observations, beta, eprob);
    calc_gammas(*model_params, alpha, beta, gamma, eprob);
    if (parameters->fixtrans == 0) {
        calc_xis(*model_params, alpha, beta, xi, eprob);
    }

    mll_diff = final_mll - initial_mll;
    bt_iterations = 0; 
    while ((mll_diff)/final_mll > TOL || (-1.0*mll_diff)/final_mll > TOL) {

        initial_mll = final_mll;
      
        /* store last */
        last_mu_ratio = (*model_params)->mu_ratio; 
        last_sig_ratio = (*model_params)->sigma_ratio; 
        last_sig_pi = (*model_params)->sigma_pi; 
        /* expected value of log(P) */
        exp_log_prob(*model_params, observations, gamma, &elogp);

        /* re-estimate */
        /* initial probabilities */
        for (i = 0; i < (*model_params)->N; i++) {
            /* (*model_params)->pi[i] = 0.0001/(*model_params)->N + 0.9999*(gamma[0][i] + 
                                            gamma[0][i + (*model_params)->N]); */
            (*model_params)->pi[i] = gamma[0][i] + gamma[0][i + (*model_params)->N];
            
        }

        /* transition probabilities and sums for emission parameters */
        total_pialt_num = 0.0;
        total_pialt_den = 0.0;
        mu_ratio_num = 0.0;
        mu_ratio_den = 0.0;

        /* first take care of transition probabilities */
        norm = 0.0;
        sum = 0.0;
        for (i = 0; i < (*model_params)->N; i++) {
            if (parameters->fixtrans == 0) {
                for (j = 0; j < (*model_params)->N; j++) {
                    for (t = 0; t < (*model_params)->T - 1; t++) {
                        addend = xi[t][i][j] + xi[t][i + (*model_params)->N][j] +
                                 xi[t][i][j + (*model_params)->N] +
                                 xi[t][i + (*model_params)->N][j + (*model_params)->N];
                        if (j != i) {
                            sum += addend;
                        }
                        norm += addend;
                    }
                }
            }

            /* state info (for estimating ratio and pi_alt parameters) */
            thisstate = &(((*model_params)->states)[i]);
            total = thisstate->total;
            minor = thisstate->minor;
            explow = (total == 0) ? 0.5 :
                      1.0*minor/total + rho_contam * (0.5 - 1.0*minor/total);
            exphigh = (total == 0) ? 0.5 :
                      1.0 - 1.0*minor/total + rho_contam * (1.0*minor/total - 0.5);
            
            for (t = 0; t < (*model_params)->T; t++) {
                mu_ratio_num += (&(observations[t]))->depthratio *
                       (gamma[t][i] + gamma[t][i + (*model_params)->N]);
                mu_ratio_den += (gamma[t][i] + gamma[t][i + (*model_params)->N]) * (rho_contam + (1.0 - rho_contam)*total*0.5);

                difflow = (&(observations[t]))->pi - explow;
                diffhigh = (&(observations[t]))->pi - exphigh;
                total_pialt_num += pow(difflow, 2) * gamma[t][i]/(explow*exphigh*4.0);
                total_pialt_num += pow(diffhigh, 2) * 
                                     gamma[t][i + (*model_params)->N] /
                                     (explow * exphigh * 4.0);
;
    
                total_pialt_den += gamma[t][i] + gamma[t][i + (*model_params)->N];
            }
        }
        if (parameters->fixtrans == 0) {
            (*model_params)->trans_prob = 0.5 * sum / (norm * (((*model_params)->N) - 1.0));
            fprintf(stderr, "New trans_prob value of %lf\n", (*model_params)->trans_prob);
            for (i = 0; i < (*model_params)->N; i++) {
                for (j = 0; j < (*model_params)->N; j++) {
                    if (i != j) {
                        ((*model_params)->a)[i][j] = (*model_params)->trans_prob;
                        ((*model_params)->a)[i+(*model_params)->N][j] = (*model_params)->trans_prob;
                        ((*model_params)->a)[i][j+(*model_params)->N] = (*model_params)->trans_prob;
                        ((*model_params)->a)[i+(*model_params)->N][j+(*model_params)->N] = (*model_params)->trans_prob;
                    }
                    else {
                        new_diag = 0.5 - (*model_params)->trans_prob * 
                                          ((*model_params)->N - 1.0);
                        ((*model_params)->a)[i][j] = new_diag;
                        ((*model_params)->a)[i+(*model_params)->N][j] = new_diag;
                        ((*model_params)->a)[i][j+(*model_params)->N] = new_diag;
                        ((*model_params)->a)[i+(*model_params)->N][j+(*model_params)->N] = new_diag;
                    }
                }
            }
        }
        (*model_params)->mu_ratio = mu_ratio_num/mu_ratio_den;
        fprintf(stderr, "New mu_ratio value of %lf\n", (*model_params)->mu_ratio);

        (*model_params)->sigma_pi = sqrt(total_pialt_num/total_pialt_den);
        fprintf(stderr, "New sigma_pi value of %lf\n", (*model_params)->sigma_pi);

        sig_ratio_num = 0.0;
        sig_ratio_den = 0.0;
        small_sig_ratio_num = 0.0;
        small_sig_ratio_den = 0.0;
        for (i = 0; i < (*model_params)->N; i++) {
            /* state info (for estimating ratio parameters) */
            thisstate = &(((*model_params)->states)[i]);
            total = thisstate->total;

            eff_cn = rho_contam + (1.0 - rho_contam) * total/2.0;
           
            for (t = 0; t < (*model_params)->T; t++) {
                addend = pow( (&(observations[t]))->depthratio -
                              eff_cn*(*model_params)->mu_ratio, 2) *
                              (gamma[t][i] + gamma[t][i + (*model_params)->N])/eff_cn;
                if (sig_ratio_num == 0.0 || addend/sig_ratio_num > 0.000000001) {
                    sig_ratio_num += addend;
                }
                else { /* aggregate small terms */
                    small_sig_ratio_num += addend;
                }

                addend = gamma[t][i] + gamma[t][i + (*model_params)->N];

                if (sig_ratio_den == 0.0 || addend/sig_ratio_den > 0.000000001) {
                    sig_ratio_den += addend;
                }
                else {
                    small_sig_ratio_den += addend;
                }
            }
        }

        /* fprintf(stderr, "Muratio num %lf + %lf, den %lf + %lf\n", sig_ratio_num, small_sig_ratio_num, sig_ratio_den, small_sig_ratio_den); */
        (*model_params)->sigma_ratio = sqrt((sig_ratio_num + small_sig_ratio_num) /
                                        (sig_ratio_den + small_sig_ratio_den));
        fprintf(stderr, "New sigma_ratio value of %lf\n", (*model_params)->sigma_ratio);

        /* This section is an optional derivative check to assure precision is still good */
        dbaum_dmuratio(*model_params, observations, gamma, &derivative);
        fprintf(stderr, "Muratio derivative is %lf\n", derivative);
        dbaum_dsigratio(*model_params, observations, gamma, &derivative);
        fprintf(stderr, "Sigratio derivative is %lf\n", derivative);
        dbaum_dsigpi(*model_params, observations, gamma, &derivative);
        fprintf(stderr, "Sigpi derivative is %lf\n", derivative);
        if (parameters->fixtrans == 0) {
            dbaum_dtransprob(*model_params, observations, xi, &derivative);
            fprintf(stderr, "Transprob derivative is %lf\n", derivative);
        }

        exp_log_prob(*model_params, observations, gamma, &elogpnew);
        if (elogpnew < elogp) {
            fprintf(stderr, "LOWER EXPECTED PROB!!!!\n");
        }
        fprintf(stderr, "New/Old expected logp values with old states: %lf, %lf\n", elogp, elogpnew);

        calc_eprobs(*model_params, observations, eprob);
        calc_alphas(*model_params, observations, alpha, eprob, &final_mll);
        calc_betas(*model_params, observations, beta, eprob);
        calc_gammas(*model_params, alpha, beta, gamma, eprob);
        if (parameters->fixtrans == 0) {
            calc_xis(*model_params, alpha, beta, xi, eprob);
        }

        mll_diff = final_mll - initial_mll;
        if (mll_diff >= 0) {
            fprintf(stderr, "MLL difference: %lf (%lf, %lf)\n", mll_diff, initial_mll, final_mll);
        }
        else {
            while (mll_diff < 0) {
                fprintf(stderr, "ERROR!!!BACKTRACKING: MLL difference: %lf (%lf, %lf)\n", mll_diff, initial_mll, final_mll);
                random_number = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
                (*model_params)->mu_ratio = last_mu_ratio + random_number;
                random_number = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
                (*model_params)->sigma_ratio = last_sig_ratio + random_number;
                random_number = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
                (*model_params)->sigma_pi = last_sig_pi + random_number;
                calc_eprobs(*model_params, observations, eprob);
                calc_alphas(*model_params, observations, alpha, eprob, &final_mll);
                calc_betas(*model_params, observations, beta, eprob);
                calc_gammas(*model_params, alpha, beta, gamma, eprob);
                if (parameters->fixtrans == 0) {
                    calc_xis(*model_params, alpha, beta, xi, eprob);
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
    free_double_matrix(eprob, (*model_params)->T);
    free_double_matrix(alpha, (*model_params)->T);
    free_double_matrix(beta, (*model_params)->T);
    free_double_matrix(gamma, (*model_params)->T);
    if (parameters->fixtrans == 0) {
        free_double_tensor(xi, (*model_params)->T, 2*(*model_params)->N);
    }

    return final_mll;
}

/* pre-calculate the probability of an observation from a particular state at a particular time */

void calc_eprobs(ModelParams *model_params, Observation *observations, double **eprob) {
    int i, t, istate, isig, minor, total;
    double logratioprob, logaltprob, sumlogs;
    double gaussnorm, rho_contam, exp_pi, exp_picontam;
    double *statenorm;
    StateData *state;

    /* extra factor for normalization of Gaussian densities */    
    gaussnorm = 1.0/(4.0*3.14159*model_params->sigma_pi*model_params->sigma_ratio);

    rho_contam = model_params->rho_contam;
    /* pre-calculate normalization factor for each different state: */
    statenorm = alloc_double_vector(model_params->N);

    for (i = 0; i < model_params->N; i++) {
        state = &((model_params->states)[i]);
        minor = state->minor;
        total = state->total;

        exp_pi = (total == 0) ? 0.5 : 1.0*minor/total;
        exp_picontam = rho_contam * (0.5 - exp_pi) + exp_pi;
        statenorm[i] = gaussnorm/sqrt(exp_picontam * (1.0 - exp_picontam));
        statenorm[i] /= sqrt(rho_contam + 0.5*(1.0 - rho_contam)*total);
        /* fprintf(stderr, "State norm for state %d is %lf\n", i, statenorm[i]); */
    }

    for (t = 0; t < model_params->T; t++) {
        for (i = 0; i < 2*model_params->N; i++) {
            istate = (i >= model_params->N) ? i - model_params->N : i;
            isig = (i >= model_params->N) ? 1 : -1;
            ratio_log_prob(model_params, observations[t], istate, &logratioprob);
            altfreq_log_prob(model_params, observations[t], istate, isig, &logaltprob);
            sumlogs = logratioprob + logaltprob;
            if ((sumlogs < 0.0000001 * abs(logratioprob) && 
                           sumlogs > -0.0000001 * abs(logratioprob)) ||
                (sumlogs < 0.0000001 * abs(logaltprob) && 
                           sumlogs > -0.0000001 * abs(logaltprob))) {
                fprintf(stderr, "Loss of precision in calculating eprobs at t=%ld for state %d/%d\n", t, istate, isig);
                exit (1); 
            }
            if (logratioprob + logaltprob < parameters->minexparg) {
                /* fprintf(stderr, "Setting eprob to 0.0 at %ld for state %d/%d (%lf, %lf)\n", t, istate, isig, logratioprob, logaltprob); */
                eprob[t][i] = 0.0;
            }
            else {
                eprob[t][i] = exp(logratioprob + logaltprob) * statenorm[istate];
            }
        }
    }
    free_double_vector(statenorm);
    /* fprintf(stderr, "Done with eprobs\n"); */
}

void calc_gammas(ModelParams *model_params, double **alpha, double **beta, double **gamma, double **eprob) {
    int i, j;
    int t;
    double norm, sum;

    for (t = 0; t < model_params->T; t++) {
        norm = 0.0;
        for (j = 0; j < 2*model_params->N; j++) {
            gamma[t][j] = alpha[t][j] * beta[t][j];
            norm += gamma[t][j];
        }
        sum = 0.0;
        for (j = 0; j < 2*model_params->N; j++) {
            gamma[t][j] /= norm;
            sum += gamma[t][j];
        }
    }
    /* fprintf(stderr, "Done with gamma, sum %lf\n", sum); */
}

void calc_xis(ModelParams *model_params, double **alpha, double **beta, double ***xi, double **eprob) {
    int i, j;
    int t;
    double norm;

    for (t = 0; t < model_params->T - 1; t++) {
        norm = 0.0;
        for (i = 0; i < 2*model_params->N; i++) {
            for (j = 0; j < 2*model_params->N; j++) {
                xi[t][i][j] = alpha[t][i] * beta[t+1][j] * model_params->a[i][j] * eprob[t+1][j];
                norm += xi[t][i][j];
            }
        }
        for (i = 0; i < 2*model_params->N; i++) {
            for (j = 0; j < 2*model_params->N; j++) {
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

    min_contam = 0.01;
    max_contam = 0.50;
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

