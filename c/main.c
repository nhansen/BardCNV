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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "bard.h"

extern Params *parameters;

void show_usage()
{
    fprintf(stderr, "Usage: bardcnv COMMAND [parameters]\n"
        "       where COMMAND is one of:\n"
        "    baumwelch - train model using observations\n"
        "    bwcontam  - train model, allowing for contamination\n"
        "    viterbi   - report most likely state at each observation\n"
        "    plotmu    - plot the probability for various values of mu (depth-ratio)\n"
        "For specific parameters for each command, run 'bardcnv COMMAND'\n");
    exit(1);
}

main(int argc, char **argv)
{
    Observation *observations;
    ModelParams *model_params;
    Observation *thisObs;
    int *states;
    double *probs, log_prob;

    get_params(argc, argv);
    small_double_check(&(parameters->minexparg));

    if (parameters->program == NULL) {
        show_usage();
    }
    if (strcmp(parameters->program, "viterbi") == 0) {
        if (!parameters->modelfile || !parameters->obsfile) {
            fprintf(stderr, "Usage: bardcnv viterbi -modelfile <file of HMM params> -obsfile <file of observation values>\n");
            exit(1);
        }
        fprintf(stderr, "Model file: %s, Observation file: %s\n", parameters->modelfile, parameters->obsfile);
        read_model(parameters->modelfile, &model_params);
        read_observations(parameters->obsfile, &model_params, &observations);
        filter_highcopy_observations(&model_params, &observations);
        states = alloc_int_vector(model_params->T);
        probs = alloc_double_vector(model_params->T);
        run_viterbi(model_params, observations, states, probs, &log_prob);
        write_states(states, probs, log_prob, model_params, observations);
        free_int_vector(states);
        free_double_vector(probs);
        fprintf(stderr, "Finished writing states.\n");
    }
    else if (strcmp(parameters->program, "baumwelch") == 0 ) {
        if (!parameters->modelfile || !parameters->obsfile) {
            fprintf(stderr, "Usage: bardcnv baumwelch -modelfile <file of HMM params> -obsfile <file of observation values>\n");
            exit(1);
        }
        fprintf(stderr, "Model file: %s, Observation file: %s\n", parameters->modelfile, parameters->obsfile);
        read_model(parameters->modelfile, &model_params);
        read_observations(parameters->obsfile, &model_params, &observations);
        filter_highcopy_observations(&model_params, &observations);
        run_baumwelch(&model_params, observations); 
        write_model(model_params);
        fprintf(stderr, "Finished writing model.\n");
    }
    else if (strcmp(parameters->program, "plotmu") == 0 ) {
        if (!parameters->modelfile || !parameters->obsfile) {
            fprintf(stderr, "Usage: bardcnv  plotmu -modelfile <file of HMM params> -obsfile <file of observation values>\n");
            exit(1);
        }
        fprintf(stderr, "Model file: %s, Observation file: %s\n", parameters->modelfile, parameters->obsfile);
        read_model(parameters->modelfile, &model_params);
        read_observations(parameters->obsfile, &model_params, &observations);
        plot_mu_prob(model_params, observations); 
    }
    else if (strcmp(parameters->program, "bwcontam") == 0) {
        if (!parameters->modelfile || !parameters->obsfile) {
            fprintf(stderr, "Usage: bardcnv bwcontam -modelfile <file of HMM params> -obsfile <file of observation values>\n");
            exit(1);
        }
        fprintf(stderr, "Model file: %s, Observation file: %s\n", parameters->modelfile, parameters->obsfile);
        read_model(parameters->modelfile, &model_params);
        read_observations(parameters->obsfile, &model_params, &observations);
        filter_highcopy_observations(&model_params, &observations);
        run_baumwelch_contam_optimization(&model_params, observations);
        write_model(model_params);
    }
    else if (strcmp(parameters->program, "precision") == 0) {
        fprintf(stdout, "Performing precision check!\n");
    }
    else {
        fprintf(stderr, "Invalid program name %s.\n", parameters->program);
        show_usage();
    }
}
