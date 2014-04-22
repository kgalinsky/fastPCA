/*
 * fastpca-iterative.c
 *
 *  Created on: Aug 8, 2013
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
size_t J = 1;
size_t K = 5;
size_t L = 10;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args(int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

int main (int argc, char **argv) {
    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");

    size_t i;
    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    kjg_geno *X  = kjg_geno_alloc(m, n);
    uint8_t *x   = malloc(sizeof(uint8_t)*n);
    double *y    = malloc(sizeof(double)*n);
    double *M    = malloc(sizeof(double)*m);

    gsl_rng *r    = kjg_rng_init();
    gsl_matrix *G1 = gsl_matrix_alloc(n, L);
    gsl_matrix *G2 = gsl_matrix_alloc(n, L);
    gsl_matrix *H  = gsl_matrix_alloc(m, (I+1)*L);
    gsl_matrix *X1;
    gsl_matrix *V1;
    gsl_vector *S1;
    gsl_vector *work1;
    gsl_matrix *Q  = gsl_matrix_alloc(m, (I+1)*L);
    gsl_matrix *T;
    gsl_matrix *X2;
    gsl_matrix *W;
    gsl_vector *S2;
    gsl_vector *work2;

    // PREP - read genotype file into memory
    kjg_genoIO_fread_geno(X, fh_geno);
    fclose(fh_geno);

    // STEP 1A - generate G - O(NL)
    kjg_matrix_set_ran_ugaussian(G1, r);
    // kjg_matrix_fprintf(fh_G, G1, "%g");
    // fclose(fh_G);

    // STEP 1B - compute H - O(MN(I+1)L)
    kjg_geno_row_means(X, M, x);
    kjg_blanczos(X, M, x, y, H, G1, G2);
    // kjg_matrix_fprintf(fh_H, H, "%g");
    // fclose(fh_H);
    gsl_matrix_free(G2);

    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    for (i = 0; i <= I; i += J) {
        fprintf(fh_evec, "#%d\n", J);
        // STEP 2P - copy H to Q
        gsl_matrix_const_view Hv = gsl_matrix_const_submatrix(H, 0, 0, H->size1, (i+1)*L);
        gsl_matrix_view Qv = gsl_matrix_submatrix(Q, 0, 0, Q->size1, (i+1)*L);
        gsl_matrix_memcpy(&Qv.matrix, &Hv.matrix);

        // STEP 2 - supposed to be pivoted QR, but can't figure it out - O(M[(I+1)L)]^2)
        X1    = gsl_matrix_alloc(Qv.matrix.size2, Qv.matrix.size2);
        V1    = gsl_matrix_alloc(Qv.matrix.size2, Qv.matrix.size2);
        S1    = gsl_vector_alloc(Qv.matrix.size2);
        work1 = gsl_vector_alloc(Qv.matrix.size2);
        gsl_linalg_SV_decomp_mod(&Qv.matrix, X1, V1, S1, work1);
        gsl_matrix_free(X1);
        gsl_matrix_free(V1);
        gsl_vector_free(S1);
        gsl_vector_free(work1);

        // STEP 3 - O(MN(I+1)L)
        T  = gsl_matrix_alloc(n, Qv.matrix.size2);
        kjg_XTH(X, M, x, y, &Qv.matrix, T);

        // STEP 4 - final SVD
        X2 = gsl_matrix_alloc(T->size2, T->size2);
        W  = gsl_matrix_alloc(T->size2, T->size2);
        S2 = gsl_vector_alloc(T->size2);
        work2 = gsl_vector_alloc(T->size2);
        gsl_linalg_SV_decomp_mod(T, X2, W, S2, work2);
        gsl_matrix_free(X2);
        gsl_vector_free(work2);
        gsl_matrix_free(W);

        gsl_matrix *V2 = T;
        gsl_matrix_view V = gsl_matrix_submatrix(V2, 0, 0, V2->size1, K);
        gsl_vector_view S = gsl_vector_subvector(S2, 0, K);
        gsl_vector_mul(&S.vector, &S.vector);

        kjg_evec_fprintf(fh_evec, &S.vector, &V.matrix, "%g");
        gsl_matrix_free(V2);
        gsl_vector_free(S2);
    }
    kjg_geno_free(X);
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
        case 'j':
            J = atoi(optarg);
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

    if (K >= L) {
        fprintf(stderr, "K >= L\n");
        exit(1);
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
