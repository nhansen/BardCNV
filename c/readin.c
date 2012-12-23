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

extern Params *parameters;

void read_observations(char *filename, ModelParams **model_params, Observation **observations)
{
    char cmd[100000];
    char nextChar;
    FILE *fc, *fp;
    char wcstring[100000];
    long noLines, obsNo;
    char chr[1000];
    char *chrptr;
    long pos;
    double ratio, pi;
    Observation *thisObs;

    /* Number of lines in observation file is upper limit on number of observations */
    sprintf(cmd, "wc -l %s | awk '{print $1}'", filename);
    if ((fc = popen(cmd, "r")) == NULL) {
        fprintf(stderr, "Error running cmd %s.\n", cmd);
        exit(1);
    }
    else { /* read number of lines and allocate memory */
        if (fscanf(fc, "%s", wcstring) == 1) {
            noLines = atoi(wcstring);
            *observations = (Observation *)malloc(sizeof(Observation)*noLines);
        }
        else {
            fprintf(stderr, "Error assessing size of file %s.\n", filename);
            exit(1);
        }
        pclose(fc);
    }

    /* Open observation file and read observations into Observation array */

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Unable to read file %s\n", filename);
        exit(1);
    }

    obsNo = 0;
    while (fscanf(fp, "%s\t%d\t%lf\t%lf", chr, &pos, &ratio, &pi) == 4) {
        chrptr = (char *)malloc(sizeof(char)*(strlen(chr) + 1));
        strcpy(chrptr, chr);
        thisObs = *observations + obsNo;
        if (ratio == 0.00) {
            pi = 0.5;
        }
        obsNo++;
        thisObs->chr = chrptr;
        thisObs->pos = pos;
        thisObs->depthratio = ratio;
        thisObs->pi = pi;
    }
    (*model_params)->T = obsNo;
    fclose(fp);
}

void read_model(char *filename, ModelParams **model_params)
{
    int noStates, i, j;
    int minor, total;
    double statepi;
    FILE *fp;
    double *piptr;
    StateData *stateptr;
    double *transProb;
    char keyword[1000];
    double offdiag, value;

    /* allocate memory for model parameters */
    *model_params = (ModelParams *)malloc(sizeof(ModelParams));
    (*model_params)->T = 0;
    /* Open observation file and read observations into Observation array */

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Unable to read file %s\n", filename);
        exit(1);
    }

    if (fscanf(fp, "N= %d", &noStates) == 1) {
        (*model_params)->N = noStates;
        (*model_params)->pi = (double *)malloc(sizeof(double)*noStates);
        (*model_params)->states = (StateData *)malloc(sizeof(StateData)*noStates);
        /* allocate space for state transitions */
        (*model_params)->a = (double **)malloc(sizeof(double *)*noStates);
        for (i=0; i<noStates; i++) {
            *((*model_params)->a + i) = (double *)malloc(sizeof(double)*noStates); 
        }
        piptr = (*model_params)->pi;
        stateptr = (*model_params)->states;
        
        for (i=0; i<noStates; i++) {
            if (fscanf(fp, "%d %d %lf", &(stateptr->minor), &(stateptr->total), piptr) == 3) {
                fprintf(stderr, "State %d %d %lf\n", stateptr->minor, stateptr->total, *piptr);
                piptr++;
                stateptr++;
            }
            else {
                fprintf(stderr, "State line must be followed by N lines of state descriptions\n");
                exit(1);
            }
        }
        fscanf(fp, " %s ", keyword);
        /* if (fscanf(fp, " Transitions ") == 0) { */
        if (strcmp(keyword, "Transitions") == 0) {
            for (i=0; i<noStates; i++) {
                transProb = (*model_params)->a[i];
                for (j=0; j<noStates; j++) {
                    if (fscanf(fp, "%lf ", transProb) == 1) {
                        transProb++;
                    }
                    else {
                        fprintf(stderr, "Transitions line must be followed by N lines of transition probabilities.\n");
                        exit(1);
                    }
                }
            }
        }
        else if (strcmp(keyword, "Transprob=") == 0) {
            if (fscanf(fp, " %lf ", &offdiag)==1) {
                fprintf(stderr, "Transprob is %lf\n", offdiag);
                for (i=0; i<noStates; i++) {
                    transProb = (*model_params)->a[i];
                    for (j=0; j<noStates; j++) {
                        if (i != j) {
                            *(transProb++) = offdiag;
                        }
                        else {
                            *(transProb++) = 1.0 - offdiag * (noStates - 1);
                        }
                    }
                }
            }
        }
        else {
            fprintf(stderr, "Transitions must follow state descriptions.\n");
            exit(1);
        }
        while (fscanf(fp, " %s %lf", keyword, &value) == 2) {
            if (strcmp(keyword, "mu_ratio=")==0) {
                (*model_params)->mu_ratio = value;
                fprintf(stderr, "mu_ratio is %lf\n", (*model_params)->mu_ratio);
            }
            else if (strcmp(keyword, "sigma_ratio=")==0) {
                (*model_params)->sigma_ratio = value;
                fprintf(stderr, "sigma_ratio is %lf\n", (*model_params)->sigma_ratio);
            }
            else if (strcmp(keyword, "rho_contam=")==0) {
                (*model_params)->rho_contam = value;
            }
            else if (strcmp(keyword, "sigma_pi=")==0) {
                (*model_params)->sigma_pi = value;
                fprintf(stderr, "sigma_pi is %lf\n", (*model_params)->sigma_pi);
            }
            else {
                fprintf(stderr, "Unrecognized keyword %s.\n", keyword);
                exit(1);
            } 
        }
    }
    else {
        fprintf(stderr, "First line of model parameters file must contain number of states!\n");
        exit(1);
    }
    fclose(fp);
}
