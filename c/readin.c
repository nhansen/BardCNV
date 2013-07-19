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
    size_t obsmem;
    char chr[1000];
    char *chrptr;
    long pos;
    long normaldepth, tumordepth, tumoralt;
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
            obsmem = sizeof(Observation)*noLines;
            *observations = (Observation *)malloc(obsmem);
            if (!*observations) {
                fprintf(stderr, "Couldn't allocate memory!\n");
                exit(1);
            }
            else {
                fprintf(stderr, "Allocated memory for observations: %lu\n", (unsigned long)obsmem);
            } 
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
    while (fscanf(fp, "%s\t%ld\t%ld\t%ld\t%ld", chr, &pos, &normaldepth, &tumordepth, &tumoralt) == 5) {
        chrptr = (char *)malloc(sizeof(char)*(strlen(chr) + 1));
        if (!chrptr) {
            fprintf(stderr, "Couldn't allocate memory!\n");
            exit(1);
        }
        strcpy(chrptr, chr);
        thisObs = *observations + obsNo;
        obsNo++;
        thisObs->chr = chrptr;
        thisObs->pos = pos;
        thisObs->normaltotaldepth = normaldepth;
        thisObs->tumortotaldepth = tumordepth;
        thisObs->tumoraltdepth = tumoralt;
    }
    (*model_params)->T = obsNo;
    fclose(fp);
}

void read_model(char *filename, ModelParams **model_params)
{
    int noStates, i, j, istate, jstate;
    int minor, total;
    double statepi;
    FILE *fp;
    StateData *stateptr;
    double *transProb;
    char keyword[1000];
    double offdiag, value;

    /* allocate memory for model parameters */
    *model_params = (ModelParams *)malloc(sizeof(ModelParams));
    if (!model_params) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    (*model_params)->T = 0;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Unable to read file %s\n", filename);
        exit(1);
    }

    if (fscanf(fp, "N= %d", &noStates) == 1) {
        (*model_params)->N = noStates;
        (*model_params)->pi = alloc_double_vector(noStates);
        (*model_params)->states = (StateData *)malloc(sizeof(StateData)*noStates);
        if (!(*model_params)->states) {
            fprintf(stderr, "Couldn't allocate memory!\n");
            exit(1);
        }
        stateptr = (*model_params)->states;
        
        for (i=0; i<noStates; i++) {
            if (fscanf(fp, "%d %d %lf", &(stateptr->minor), &(stateptr->total), &(*((*model_params)->pi + i))) == 3) {
                fprintf(stderr, "State %d %d %lf\n", stateptr->minor, stateptr->total, (*model_params)->pi[i]);
                stateptr++;
            }
            else {
                fprintf(stderr, "State line must be followed by N lines of state descriptions\n");
                exit(1);
            }
        }
        /* allocate space for state transitions */
        (*model_params)->a = alloc_double_matrix(noStates, noStates);
        fscanf(fp, " %s ", keyword);
        /* if (strcmp(keyword, "Transitions") == 0) {
            for (i=0; i < noStates; i++) {
                transProb = (*model_params)->a[i];
                for (j=0; j < noStates; j++) {
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
        else */
        if (strcmp(keyword, "Transprob=") == 0) {
            if (fscanf(fp, " %lf ", &((*model_params)->trans_prob))==1) {
                fprintf(stderr, "Transprob is %g\n", (*model_params)->trans_prob);
                for (i=0; i < noStates; i++) {
                    transProb = (*model_params)->a[i];
                    for (j=0; j < noStates; j++) {
                        istate = (i >= noStates) ? i - noStates : i;
                        jstate = (j >= noStates) ? j - noStates : j;
                        if (istate != jstate) {
                            *(transProb++) = (*model_params)->trans_prob;
                        }
                        else {
                            *(transProb++) = 1.0 - (*model_params)->trans_prob * (noStates - 1.0);
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
            if (strcmp(keyword, "mu=")==0) {
                (*model_params)->mu = value;
                fprintf(stderr, "mu is %lf\n", (*model_params)->mu);
            }
            else if (strcmp(keyword, "sigma=")==0) {
                (*model_params)->sigma = value;
                fprintf(stderr, "sigma is %lf\n", (*model_params)->sigma);
            }
            else if (strcmp(keyword, "rho_contam=")==0) {
                (*model_params)->rho_contam = value;
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

void filter_highcopy_observations(ModelParams **model_params, Observation **observations)
{
    int i, j, maxcopies, totalcopies, noStates, noskipped;
    double maxdepth, thistumordepth, thisnormaldepth;

    if (parameters->maxratio == 0) {
        noStates = (*model_params)->N;
        maxcopies = 0;
        for (i=0; i < noStates; i++) {
            totalcopies = (&(((*model_params)->states)[i]))->total;
            if (totalcopies > maxcopies) {
                maxcopies = totalcopies;
            }
        }
        maxdepth = maxcopies*(*model_params)->mu*0.5 + 3.0*(*model_params)->sigma;
    }
    else {
        maxdepth = parameters->maxratio;
    }

    fprintf(stderr, "Max depth is %lf\n", maxdepth);
    noskipped = 0;
    j = 0; /* index of observations slot to be filled */
    for (i=0; i < (*model_params)->T; i++) {
        thistumordepth = (&(*observations)[i])->tumortotaldepth;
        thisnormaldepth = (&(*observations)[i])->normaltotaldepth;
        if (thistumordepth > maxdepth*thisnormaldepth) {
            /* fprintf(stderr, "Skipping observation with depth %lf\n", thisdepth); */
            noskipped++;
            continue;
        }
        else {
            (*observations)[j] = (*observations)[i]; 
            j++;
        }
    }
    (*model_params)->T -= noskipped;
    if (noskipped > 0) {
        fprintf(stderr, "Skipped %i observations with coverage ratios higher than %lf.\n", noskipped, maxdepth);
    }
}
