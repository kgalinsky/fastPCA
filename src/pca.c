/*
 * pca.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <lapacke.h>

#include "kjg_geno.h"
#include "kjg_genoIO.h"
#include "kjg_GRM.h"
#include "kjg_util.h"
#include "kjg_gsl.h"

// options and arguments
size_t K = 10;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args (int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

// helper functions
void scale_evals (gsl_vector* evals);

// MAIN FUNCTION
int main (int argc, char** argv) {
    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");
    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    kjg_geno* X = kjg_geno_alloc(m, n);
    double* M = malloc(m * sizeof(double));

    kjg_genoIO_fread_geno(X, fh_geno);
    kjg_geno_row_means(X, M);

    kjg_GRM* GRM = kjg_GRM_alloc(X);
    kjg_GRM_calc(GRM, X, M);

    gsl_vector* evals = gsl_vector_alloc(n);
    gsl_matrix* evecs = gsl_matrix_alloc(n, n);

    LAPACKE_dspevd(LAPACK_COL_MAJOR, 'V', 'L', n, GRM->data, evals->data,
            evecs->data, n);

    kjg_GRM_free(GRM);

    scale_evals(evals);
    gsl_vector_reverse(evals);

    {
        size_t i;

        for (i = 0; i < n/2; i++)
            gsl_matrix_swap_columns(evecs, i, n-i-1);

    }

    gsl_vector_view evalsk = gsl_vector_subvector(evals, 0, K);
    gsl_matrix_view evecsk = gsl_matrix_submatrix(evecs, 0, 0, n, K);

    kjg_gsl_evec_fprintf(fh_evec, &evalsk.vector, &evecsk.matrix, "%g");

    gsl_vector_free(evals);
    gsl_matrix_free(evecs);

    return (0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:")) != -1) {
        switch (c) {
        case 'k':
            K = atoi(optarg);
            break;
        case 'o':
            OUTPUT_PREFIX = optarg;
            break;
        default:
            fprintf(stderr, "Unrecognized option '-%c'\n", optopt);
            exit(1);
        }
    }

    if (optind == argc - 1) {
        GENO_FILENAME = argv[optind];
    }
    else if (optind == argc) {
        fprintf(stderr, "No input file specified\n");
        exit(1);
    }
    else {
        fprintf(stderr, "Too many arguments\n");
        exit(1);
    }

    if (OUTPUT_PREFIX == NULL) {
        OUTPUT_PREFIX = GENO_FILENAME;
    }
}

void scale_evals (gsl_vector* evals) {
    double sum = 0;
    double *p = evals->data;
    size_t i = evals->size;

    while (i--)
        sum += *(p++);
    double scale = (double) evals->size / sum;

    gsl_vector_scale(evals, scale);
}
