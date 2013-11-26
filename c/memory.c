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

double **alloc_double_matrix(int rows, int cols) {

    int i;
    double **matrix;

    /* this routine will allocate rows*cols doubles and rows pointers to doubles */
    matrix = (double **)malloc(sizeof(double *) * rows); 
    if (!matrix) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(sizeof(double) * cols);
        if (!matrix[i]) {
            fprintf(stderr, "Couldn't allocate memory!\n");
            exit(1);
        }
    }

    return matrix;
}

void free_double_matrix(double **matrix, int rows) {

    int i;

    for (i = 0; i < rows; i++) {
        free(matrix[i]);
    }

    free(matrix);
}

double ***alloc_double_tensor(int rows, int cols, int levels) {

    int i, j;
    double ***tensor;

    /* this routine will allocate rows*cols*levels doubles */
    tensor = (double ***)malloc(sizeof(double **) * rows); 
    if (!tensor) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        *(tensor + i) = (double **)malloc(sizeof(double *) * cols);
        if (!tensor[i]) {
            fprintf(stderr, "Couldn't allocate memory!\n");
            exit(1);
        }
        for (j = 0; j < cols; j++) {
            (*(tensor + i))[j] = (double *)malloc(sizeof(double) * levels);
            if (!tensor[i][j]) {
                fprintf(stderr, "Couldn't allocate memory!\n");
                exit(1);
            }
        }
    }

    return tensor;
}

void free_double_tensor(double ***tensor, int rows, int cols) {

    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            free(tensor[i][j]);
        }
        free(tensor[i]);
    }

    free(tensor);

}

double *alloc_double_vector(int entries) {

    double *vector;

    vector = (double *)malloc(sizeof(double) * entries);
    if (!vector) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    return vector;
}

void free_double_vector(double *vector) {

    free(vector);
}

int *alloc_int_vector(int entries) {

    int *vector;

    vector = (int *)malloc(sizeof(int) * entries);
    if (!vector) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    return vector;
}

void free_int_vector(int *vector) {

    free(vector);
}

int **alloc_int_matrix(int rows, int cols) {

    int i;
    int **matrix;

    matrix = (int **)malloc(sizeof(int *) * rows); 
    if (!matrix) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        *(matrix + i) = (int *)malloc(sizeof(int) * cols);
        if (!matrix[i]) {
            fprintf(stderr, "Couldn't allocate memory!\n");
            exit(1);
        }
    }

    return matrix;
}

void free_int_matrix(int **matrix, int rows) {

    int i;

    for (i = 0; i < rows; i++) {
        free(*(matrix + i));
    }

    free(matrix);
}

void copy_params(ModelParams *from_model_params, ModelParams **to_model_params) {
    int i, j, noStates, newminor, newtotal;
    StateData *stateptr;

    /* allocate memory for model parameters */
    *to_model_params = (ModelParams *)malloc(sizeof(ModelParams));
    if (!to_model_params) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }
    (*to_model_params)->T = from_model_params->T;
    (*to_model_params)->N = from_model_params->N;
    noStates = from_model_params->N;
    (*to_model_params)->trans_prob = from_model_params->trans_prob;
    (*to_model_params)->mu = from_model_params->mu;
    (*to_model_params)->sigma = from_model_params->sigma;
    (*to_model_params)->rho_contam = from_model_params->rho_contam;

    /* allocate space for state transitions, state parameters */
    (*to_model_params)->pi = alloc_double_vector(noStates);
    for (i=0; i<noStates; i++) {
        /* fprintf(stderr, "Copying from %lf\n", *(from_model_params->pi + i)); */
        *((*to_model_params)->pi + i) = *(from_model_params->pi + i);
    }
    /* fprintf(stderr, "Populated pi vector\n"); */

    /* fprintf(stderr, "Allocating %d by %d matrix\n", 2*noStates, 2*noStates); */
    (*to_model_params)->a = alloc_double_matrix(2*noStates, 2*noStates);

    for (i = 0; i < 2*noStates; i++) {
        /* transitions */
        for (j = 0; j < 2*noStates; j++) {
            (*to_model_params)->a[i][j] = from_model_params->a[i][j];
        }
    }
    (*to_model_params)->states = (StateData *)malloc(sizeof(StateData)*noStates);
    if (!(*to_model_params)->states) {
        fprintf(stderr, "Couldn't allocate memory!\n");
        exit(1);
    }

    stateptr = (*to_model_params)->states;
        
    for (i=0; i<noStates; i++) {
        stateptr->minor = (&(from_model_params->states[i]))->minor;
        stateptr->total = (&(from_model_params->states[i]))->total;
        stateptr++;
    }
    /* check */
    for (i=0; i<noStates; i++) {
        newminor = (&((*to_model_params)->states[i]))->minor;
        if (newminor != (&(from_model_params->states[i]))->minor) {
            fprintf(stderr, "Unequal minor\n");
        }
        newtotal = (&((*to_model_params)->states[i]))->total;
        if (newtotal != (&(from_model_params->states[i]))->total) {
            fprintf(stderr, "Unequal total\n");
        }
    }
}
