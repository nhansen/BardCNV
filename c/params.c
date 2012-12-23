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

Params *parameters;

void get_params(argc, argv)
    int argc;
    char **argv;
{
    int i;

    char empty[10] = "";
    parameters = (Params *) malloc(sizeof(Params));

    /* defaults */
    parameters->min = 0.0;
    parameters->max = 0.0;
    parameters->nobins = 100;

    if (argc >= 1) {
        parameters->program = (char *) malloc((strlen(argv[1]) + 1) * sizeof(char));
        strcpy(parameters->program, argv[1]);
    }
    else {
        parameters->program = empty;
    }

    for (i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-modelfile")) {
            parameters->modelfile = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->modelfile, argv[i]);
        }
        if (!strcmp(argv[i], "-obsfile")) {
            parameters->obsfile = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->obsfile, argv[i]);
        }
        if (!strcmp(argv[i], "-min")) {
            parameters->min = atof(argv[++i]);
            fprintf(stderr, "Parsed min value %lf\n", parameters->min);
        }
        if (!strcmp(argv[i], "-max")) {
            parameters->max = atof(argv[++i]);
            fprintf(stderr, "Parsed max value %lf\n", parameters->max);
        }
        if (!strcmp(argv[i], "-nobins")) {
            parameters->nobins = atoi(argv[++i]);
        }
    }
}

void copy_params(ModelParams *from_model_params, ModelParams **to_model_params) {
    int i, j, noStates;
    double *piptr, **transProb;
    StateData *stateptr;

    /* allocate memory for model parameters */
    *to_model_params = (ModelParams *)malloc(sizeof(ModelParams));
    (*to_model_params)->T = from_model_params->T;
    noStates = from_model_params->N;
    (*to_model_params)->N = noStates;
    (*to_model_params)->pi = (double *)malloc(sizeof(double)*noStates);
    (*to_model_params)->states = (StateData *)malloc(sizeof(StateData)*noStates);

    /* allocate space for state transitions */
    (*to_model_params)->a = (double **)malloc(sizeof(double *)*noStates);
    for (i=0; i<noStates; i++) {
        *((*to_model_params)->a + i) = (double *)malloc(sizeof(double)*noStates); 
    }
    piptr = (*to_model_params)->pi;
    stateptr = (*to_model_params)->states;
        
    for (i=0; i<noStates; i++) {
        *piptr = (from_model_params->pi)[i];
        piptr++;
        stateptr->minor = (&(from_model_params->states)[i])->minor;
        stateptr->total = (&(from_model_params->states)[i])->total;
        stateptr++;

        /* transitions */
        transProb = (*to_model_params)->a;
        for (j=0; j<noStates; j++) {
            transProb[i][j] = (from_model_params->a)[i][j];
        }
    }
    (*to_model_params)->mu_ratio = from_model_params->mu_ratio;
    (*to_model_params)->sigma_ratio = from_model_params->sigma_ratio;
    (*to_model_params)->sigma_pi = from_model_params->sigma_pi;
    (*to_model_params)->rho_contam = from_model_params->rho_contam;
}
