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

    gsl_rng *r    = kjg_rng_init();
    gsl_matrix *G  = gsl_matrix_alloc(n, L);
    gsl_matrix *H  = gsl_matrix_alloc(m, (I+1)*L);
    gsl_matrix *X1 = gsl_matrix_alloc(H->size2, H->size2);
    gsl_matrix *V1 = gsl_matrix_alloc(H->size2, H->size2);
    gsl_vector *S1 = gsl_vector_alloc(H->size2);
    gsl_vector *work1 = gsl_vector_alloc(H->size2);
    gsl_matrix *Q  = H; // for clarity
    gsl_matrix *T  = gsl_matrix_alloc(n, H->size2);
    gsl_matrix *X2 = gsl_matrix_alloc(T->size2, T->size2);
    gsl_matrix *W  = gsl_matrix_alloc(T->size2, T->size2);
    gsl_vector *S2 = gsl_vector_alloc(T->size2);
    gsl_vector *work2 = gsl_vector_alloc(T->size2);
    gsl_matrix *V2 = T; // for clarity

    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    // PREP - read genotype file into memory
    kjg_genoIO_fread_geno(X, fh_geno);
    fclose(fh_geno);

    // STEP 1A - generate G - O(NL)
    kjg_matrix_set_ran_ugaussian(G, r);
    // kjg_matrix_fprintf(fh_G, G1, "%g");
    // fclose(fh_G);

    // STEP 1B - compute H - O(MN(I+1)L)
    kjg_geno_row_means(X, M);
    kjg_blanczos(X, M, G, H);
    // kjg_matrix_fprintf(fh_H, H, "%g");
    // fclose(fh_H);

    // STEP 2 - supposed to be pivoted QR, but can't figure it out - O(M[(I+1)L)]^2)
    gsl_linalg_SV_decomp_mod(H, X1, V1, S1, work1);
    gsl_matrix_free(X1);
    gsl_matrix_free(V1);
    gsl_vector_free(S1);
    gsl_vector_free(work1);

    // STEP 3 - O(MN(I+1)L)
    kjg_XTH(X, M, Q, T);
    kjg_geno_free(X);

    // STEP 4 - final SVD
    gsl_linalg_SV_decomp_mod(T, X2, W, S2, work2);
    gsl_matrix_free(X2);
    gsl_vector_free(work2);
    gsl_matrix_free(W);

    gsl_matrix_view V = gsl_matrix_submatrix(V2, 0, 0, V2->size1, K);
    gsl_vector_view S = gsl_vector_subvector(S2, 0, K);
    gsl_vector_mul(&S.vector, &S.vector);
    gsl_vector_scale(&S.vector, 1.0 / m);

    kjg_evec_fprintf(fh_evec, &S.vector, &V.matrix, "%g");
    gsl_matrix_free(V2);
    gsl_vector_free(S2);
    fclose(fh_evec);

    return(0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:")) != -1) {
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
