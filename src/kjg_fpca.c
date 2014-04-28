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

void kjg_fpca_blanczos(
        const kjg_geno* X, const double* M, gsl_matrix* G,
        gsl_matrix* H) {
    size_t i;

    gsl_matrix* G2 = gsl_matrix_alloc(G->size1, G->size2);
    gsl_matrix* Gswap;
    gsl_matrix_view Hsub;

    for (i = 0; i < H->size2; i += G->size2) {
        Hsub = gsl_matrix_submatrix(H, 0, i, H->size1, G->size2);
        kjg_fpca_XTXG(X, M, G, &Hsub.matrix, G2);
        kjg_gsl_matrix_frobenius_normalize(&Hsub.matrix);

        Gswap = G2;
        G2 = G;
        G = Gswap;
    }

    gsl_matrix_free(G2);
}

void kjg_fpca_XTXG(
        const kjg_geno *X, const double* M, const gsl_matrix *G1,
        gsl_matrix *H, gsl_matrix *G2) {
    size_t i;                                   // row index
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n); // genotypes
    double *y   = malloc(sizeof(double)*X->n);  // normalized

    gsl_vector_view Hrow;
    gsl_vector_view Xrow = gsl_vector_view_array(y, X->n);

    gsl_matrix_set_zero(H);
    gsl_matrix_set_zero(G2);

    for (i = 0; i < X->m; i++) {
        kjg_geno_get_row(x, X, i);
        kjg_geno_normalize_m(M[i], x, y, X->n);
        Hrow = gsl_matrix_row(H, i);
        gsl_blas_dgemv(CblasTrans, 1, G1, &Xrow.vector, 0, &Hrow.vector);
        gsl_blas_dger (1, &Xrow.vector, &Hrow.vector, G2);
    }

    free(x);
    free(y);
}

void kjg_fpca_XG(const kjg_geno *X, const double *M, const gsl_matrix *G,
        gsl_matrix *H) {
    size_t i;                                   // row index
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n); // genotypes
    double *y   = malloc(sizeof(double)*X->n);  // normalized

    gsl_vector_view Hrow;
    gsl_vector_view Xrow = gsl_vector_view_array(y, X->n);

    gsl_matrix_set_zero(H);

    for (i = 0; i < X->m; i++) {
        kjg_geno_get_row(x, X, i);
        kjg_geno_normalize_m(M[i], x, y, X->n);
        Hrow = gsl_matrix_row(H, i);
        gsl_blas_dgemv(CblasTrans, 1, G, &Xrow.vector, 0, &Hrow.vector);
    }

    free(x);
    free(y);
}

void kjg_fpca_XTH(const kjg_geno *X, const double *M, const gsl_matrix *H,
        gsl_matrix *G) {
    size_t i;                                   // row index
    uint8_t *x  = malloc(sizeof(uint8_t)*X->n); // genotypes
    double *y   = malloc(sizeof(double)*X->n);  // normalized

    gsl_vector_view Xrow = gsl_vector_view_array(y, X->n);

    gsl_matrix_set_zero(G);

    for (i = 0; i < X->m; i++) {
        kjg_geno_get_row(x, X, i);
        kjg_geno_normalize_m(M[i], x, y, X->n);
        gsl_vector_const_view Hrow = gsl_matrix_const_row(H, i);
        gsl_blas_dger (1, &Xrow.vector, &Hrow.vector, G);
    }

    free(x);
    free(y);
}
