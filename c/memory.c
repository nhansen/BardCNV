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

    matrix = (double **)malloc(sizeof(double *) * rows); 
    for (i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(sizeof(double) * cols);
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

    tensor = (double ***)malloc(sizeof(double **) * rows); 
    for (i = 0; i < rows; i++) {
        *(tensor + i) = (double **)malloc(sizeof(double *) * cols);
        for (j = 0; j < cols; j++) {
            (*(tensor + i))[j] = (double *)malloc(sizeof(double) * levels);
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
}

void free_double_vector(double *vector) {

    free(vector);
}

int **alloc_int_matrix(int rows, int cols) {

    int i;
    int **matrix;

    matrix = (int **)malloc(sizeof(int *) * rows); 
    for (i = 0; i < rows; i++) {
        *(matrix + i) = (int *)malloc(sizeof(int) * cols);
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

