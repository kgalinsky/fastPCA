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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#include "kjg_geno.h"
#include "kjg_genoIO.h"
#include "kjg_gsl.h"
#include "kjg_util.h"

// options and arguments
size_t K = 10;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args(int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

int main(int argc, char** argv) {
    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");
    FILE *fh_eval = kjg_fopen_suffix(OUTPUT_PREFIX, "eval", "w");
    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    size_t i;

    // allocate everything
    char* buffer = malloc(sizeof(char)*(n+1));
    uint8_t* x   = malloc(sizeof(uint8_t)*n);
    double* C   = kjg_geno_correlation_matrix_init(n);
    gsl_matrix_view A = gsl_matrix_view_array(C, n, n);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    gsl_vector_view eval_K = gsl_vector_subvector(eval, 0, K);
    gsl_matrix_view evec_K = gsl_matrix_submatrix(evec, 0, 0, n, K);

    for (i = 0; i < m; i++) {
        kjg_genoIO_fread(buffer, x, n, fh_geno);
        kjg_geno_snp_correlation(x, C, n);
    }
    kjg_geno_correlation_matrix_finish(C, n);

    free(buffer);
    free(x);

    // kjg_matrix_fprintf(fh_corr, &A.matrix, "%g");

    gsl_eigen_symmv(&A.matrix, eval, evec, w);

    free(C);
    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    {
        double sum;
        for (i = 0; i < n; i++) {
            sum += gsl_vector_get(eval, i);
        }
        gsl_vector_scale(eval, (double) n / sum );
    }

    gsl_vector_fprintf(fh_eval, eval, "%g");
    kjg_gsl_evec_fprintf(fh_evec, &eval_K.vector, &evec_K.matrix, "%g");

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    fclose(fh_evec);

    return(0);
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

