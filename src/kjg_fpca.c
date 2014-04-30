/*
 * kjg_fpca.c
 *
 *  Created on: Apr 28, 2014
 *      Author: Kevin
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

#include "kjg_fpca.h"
#include "kjg_geno.h"
#include "kjg_gsl.h"

size_t KJG_FPCA_ROWS = 256;

void kjg_fpca_blanczos(
        const kjg_geno* X, const double* M, gsl_matrix* G,
        gsl_matrix* H) {
    size_t i;

    gsl_matrix* G2 = gsl_matrix_alloc(G->size1, G->size2);
    gsl_matrix* Gswap;
    gsl_matrix_view Hsub;

    for (i = 0; i < H->size2 - G->size2; i += G->size2) {
        Hsub = gsl_matrix_submatrix(H, 0, i, H->size1, G->size2);
        kjg_fpca_XTXG(X, M, G, &Hsub.matrix, G2);
        kjg_gsl_matrix_frobenius_normalize(&Hsub.matrix);

        Gswap = G2;
        G2 = G;
        G = Gswap;
    }

    gsl_matrix_free(G2);

    Hsub = gsl_matrix_submatrix(H, 0, i, H->size1, G->size2);
    kjg_fpca_XG(X, M, G, &Hsub.matrix);
    kjg_gsl_matrix_frobenius_normalize(&Hsub.matrix);
}

void kjg_fpca_XTXG(
        const kjg_geno *X, const double* M, const gsl_matrix *G1,
        gsl_matrix *H, gsl_matrix *G2) {
    size_t i, r;                                                // row index
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n);                 // genotypes
    double *Y   = malloc(sizeof(double)*X->n*KJG_FPCA_ROWS);    // normalized
    gsl_matrix_view Hmat, Xmat;

    gsl_matrix_set_zero(H);
    gsl_matrix_set_zero(G2);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        Hmat = gsl_matrix_submatrix (H, i, 0, r, H->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xmat.matrix, G1, 0, &Hmat.matrix);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xmat.matrix, &Hmat.matrix, 1, G2);
    }

    free(x);
    free(Y);
}

void kjg_fpca_XG(const kjg_geno *X, const double *M, const gsl_matrix *G,
        gsl_matrix *H) {
    size_t i, r;
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n);
    double *Y   = malloc(sizeof(double)*X->n*KJG_FPCA_ROWS);
    gsl_matrix_view Hmat, Xmat;

    gsl_matrix_set_zero(H);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        Hmat = gsl_matrix_submatrix (H, i, 0, r, H->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xmat.matrix, G, 0, &Hmat.matrix);
    }

    free(x);
    free(Y);
}

void kjg_fpca_XTH(const kjg_geno *X, const double *M, const gsl_matrix *H,
        gsl_matrix *G) {
    size_t i, r;
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n);
    double *Y   = malloc(sizeof(double)*X->n*KJG_FPCA_ROWS);
    gsl_matrix_view Xmat;

    gsl_matrix_set_zero(G);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        gsl_matrix_const_view Hmat   = gsl_matrix_const_submatrix (H, i, 0, r, H->size2);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xmat.matrix, &Hmat.matrix, 1, G);
    }

    free(x);
    free(Y);
}
