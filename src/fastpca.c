/*
 * fastpca.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_linalg.h"

#include "kjg_geno.h"
#include "kjg_genoIO.h"
#include "kjg_gsl.h"
#include "kjg_util.h"
#include "kjg_fpca.h"

//#include <lapacke.h>

// options and arguments
size_t I = 10;
size_t L = 10;
size_t K = 5;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args(int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

int main (int argc, char **argv) {
    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    kjg_geno *X  = kjg_geno_alloc(m, n);
    double *M    = malloc(sizeof(double)*m);

    gsl_rng *r    = kjg_gsl_rng_init();
    gsl_matrix *G  = gsl_matrix_alloc(n, L);
    gsl_matrix *H  = gsl_matrix_alloc(m, (I+1)*L);

    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    // PREP - read genotype file into memory
    kjg_genoIO_fread_geno(X, fh_geno);
    fclose(fh_geno);

    // STEP 1A - generate G - O(NL)
    kjg_gsl_matrix_set_ran_ugaussian(G, r);

    // STEP 1B - compute H - O(MN(I+1)L)
    kjg_geno_row_means(X, M);
    kjg_fpca_blanczos(X, M, G, H);

    // STEP 2 - supposed to be pivoted QR, but can't figure it out - O(M[(I+1)L)]^2)
    {
        size_t lwork = 5*(H->size1 + H->size2);
        double *S    = malloc(sizeof(double)*H->size2);
        double *work = malloc(sizeof(double)*lwork);
        double U, V;
        int info = LAPACKE_dgesvd(101, 'O', 'N',
                H->size1, H->size2, H->data, H->tda,
                S, &U, m, &V, m, work);
        free(S);
        free(work);
    }

    // STEP 3 - O(MN(I+1)L)
    gsl_matrix *Q  = H; // for clarity
    gsl_matrix *T  = gsl_matrix_alloc(n, H->size2);
    kjg_fpca_XTH(X, M, Q, T);
    kjg_geno_free(X);
    free(M);

    // STEP 4 - final SVD
    gsl_vector *S = gsl_vector_alloc(T->size2);
    {
        size_t lwork = 5*(T->size1 + T->size2);
        double *work = malloc(sizeof(double)*lwork);
        double U, V;
        int info = LAPACKE_dgesvd(101, 'O', 'N',
                T->size1, T->size2, T->data, T->tda, S->data,
                &U, m, &V, m, work);
        free(work);
    }

    // STEP 5 - output top K principle components
    gsl_matrix_view Vk = gsl_matrix_submatrix(T, 0, 0, T->size1, K);
    gsl_vector_view Sk = gsl_vector_subvector(S, 0, K);
    gsl_vector_mul(&Sk.vector, &Sk.vector);
    gsl_vector_scale(&Sk.vector, 1.0 / m);

    kjg_gsl_evec_fprintf(fh_evec, &Sk.vector, &Vk.matrix, "%g");
    gsl_matrix_free(T);
    gsl_vector_free(S);
    fclose(fh_evec);

    return(0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:r:")) != -1) {
        switch (c) {
        case 'i':
            I = atoi(optarg);
            break;
        case 'k':
            K = atoi(optarg);
            break;
        case 'l':
            L = atoi(optarg);
            break;
        case 'o':
            OUTPUT_PREFIX = optarg;
            break;
        case 'r':
            KJG_FPCA_ROWS = atoi(optarg);
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

    if (K >= L) {
        fprintf(stderr, "K >= L\n");
        exit(1);
    }

    if (OUTPUT_PREFIX == NULL) {
        OUTPUT_PREFIX = GENO_FILENAME;
    }
}
