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
    char *emptyptr;

    parameters = (Params *) malloc(sizeof(Params));

    emptyptr = (char *) malloc(sizeof(char));
    strcpy(emptyptr, "");

    /* defaults */
    parameters->min = 0.01; /* for plotting log likelihoods */
    parameters->max = 0.5; /* for plotting log likelihoods */
    parameters->nobins = 100; /* for plotting log likelihoods */
    parameters->fixtrans = 0; /* option to skip optimization of transition probs */
    parameters->derivatives = 0; /* option to calculate derivatives with respect to parameters */
    parameters->maxratio = 0; /* option to limit observations to those having a given maximum ratio of tumortotaldepth to normaltotaldepth */
    parameters->verbose = 0; /* option to print more output for diagnosis */
    parameters->success = 0.0; /* option to pass a success probability to "testbin" */
    parameters->region = emptyptr;
    parameters->bedfile = emptyptr;
    parameters->minqual = 0;
    parameters->mapqual = 0;

    if (argc > 1 && argv[1][0] != '-') {
        parameters->program = (char *) malloc((strlen(argv[1]) + 1) * sizeof(char));
        strcpy(parameters->program, argv[1]);
    }
    else {
        parameters->program = NULL;
    }

    for (i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-modelfile")) {
            parameters->modelfile = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->modelfile, argv[i]);
        }
        else if (!strcmp(argv[i], "-obsfile")) {
            parameters->obsfile = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->obsfile, argv[i]);
        }
        else if (!strcmp(argv[i], "-min")) {
            parameters->min = atof(argv[++i]);
            fprintf(stderr, "Parsed Min value %lf\n", parameters->min);
            if (parameters->min < 0.01) {
                parameters->min = 0;
            }
        }
        else if (!strcmp(argv[i], "-max")) {
            parameters->max = atof(argv[++i]);
            fprintf(stderr, "Parsed Max value %lf\n", parameters->max);
        }
        else if (!strcmp(argv[i], "-maxratio")) {
            parameters->maxratio = atof(argv[++i]);
            fprintf(stderr, "Parsed maxratio value %lf\n", parameters->maxratio);
        }
        else if (!strcmp(argv[i], "-nobins")) {
            parameters->nobins = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-fixtrans")) {
            parameters->fixtrans = 1;
        }
        else if (!strcmp(argv[i], "-derivatives")) {
            parameters->derivatives = 1;
        }
        else if (!strcmp(argv[i], "-verbose")) {
            parameters->verbose = 1;
        }
        else if (!strcmp(argv[i], "-success")) {
            parameters->success = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-fasta")) {
            parameters->fasta = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char));
            strcpy(parameters->fasta, argv[i]);
        }
        else if (!strcmp(argv[i], "-bam")) {
            parameters->bam = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->bam, argv[i]);
        }
        else if (!strcmp(argv[i], "-region")) {
            parameters->region = (char *)malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->region, argv[i]);
        }
        else if (!strcmp(argv[i], "-minqual")) {
            parameters->minqual = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-mapqual")) {
            parameters->mapqual = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-bedfile")) {
            parameters->bedfile = (char *)malloc((strlen(argv[++i]) + 1) * sizeof(char));
            strcpy(parameters->bedfile, argv[i]);
        }
    }
}
