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

void kjg_fpca (
        const kjg_geno* X,
        const double* M,
        gsl_vector* eval,
        gsl_matrix* evec,
        size_t L,
        size_t I) {

    if (evec->size1 != X->n) exit(1);
    if (eval->size != evec->size2) exit(1);
    if (eval->size >= L) exit(1);
    if (I == 0) exit(1);

    // PART A - compute Q such that X ~ Q * (Q^T) * X
    gsl_matrix* G1 = gsl_matrix_alloc(X->n, L);
    gsl_matrix* G2 = gsl_matrix_alloc(X->n, L);
    gsl_matrix* Q = gsl_matrix_alloc(X->m, (I + 1) * L);
    gsl_matrix* Gswap;

    gsl_rng *r = kjg_gsl_rng_init();
    kjg_gsl_matrix_set_ran_ugaussian(G1, r);

    size_t i;
    for (i = 0; i < I; i++) {
        gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, i * L, X->m, L);

        // do the multiplication
        kjg_fpca_XTXA(X, M, G1, &Qi.matrix, G2);

        // orthonormalize (Gram-Schmidt equivalent)
        kjg_gsl_matrix_QR(G2);

        Gswap = G2;
        G2 = G1;
        G1 = Gswap;
    }

    gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, I * L, X->m, L);
    kjg_fpca_XA(X, M, G1, &Qi.matrix);
    kjg_gsl_matrix_QR(Q);

    gsl_matrix_free(G1);
    gsl_matrix_free(G2);

    // PART B - compute B matrix, take SVD and return
    gsl_matrix* B = gsl_matrix_alloc(X->n, (I + 1) * L);
    kjg_fpca_XTB(X, M, Q, B);

    gsl_matrix* Vtilda = gsl_matrix_alloc(X->n, (I + 1) * L);
    gsl_matrix* Utilda = gsl_matrix_alloc((I + 1) * L, (I + 1) * L);
    gsl_vector* Stilda = gsl_vector_alloc((I + 1) * L);

    kjg_gsl_matrix_SVD(B, Utilda, Vtilda, Stilda);

    gsl_matrix_view Vk = gsl_matrix_submatrix(Vtilda, 0, 0, X->n, eval->size);
    gsl_matrix_memcpy(evec, &Vk.matrix);

    gsl_vector_view Sk = gsl_vector_subvector(Stilda, 0, eval->size);
    gsl_vector_mul(&Sk.vector, &Sk.vector);
    gsl_vector_scale(&Sk.vector, 1.0 / X->m);
    gsl_vector_memcpy(eval, &Sk.vector);

    gsl_matrix_free(Q);
    gsl_matrix_free(B);
    gsl_matrix_free(Utilda);
    gsl_matrix_free(Vtilda);
    gsl_vector_free(Stilda);
}

void kjg_svd (
        const kjg_geno* X,
        const double* M,
        gsl_vector* S,
        gsl_matrix* U,
        gsl_matrix* V,
        size_t L,
        size_t I) {

    if (U->size1 != X->m) exit(1);
    if (V->size1 != X->n) exit(1);
    if (S->size != U->size2) exit(1);
    if (S->size != V->size2) exit(1);
    if (S->size >= L) exit(1);
    if (I == 0) exit(1);

    // PART A - compute Q such that X ~ Q * (Q^T) * X
    gsl_matrix* G1 = gsl_matrix_alloc(X->n, L);
    gsl_matrix* G2 = gsl_matrix_alloc(X->n, L);
    gsl_matrix* Q = gsl_matrix_alloc(X->m, (I + 1) * L);
    gsl_matrix* Gswap;

    gsl_rng *r = kjg_gsl_rng_init();
    kjg_gsl_matrix_set_ran_ugaussian(G1, r);

    size_t i;
    for (i = 0; i < I; i++) {
        gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, i * L, X->m, L);

        // do the multiplication
        kjg_fpca_XTXA(X, M, G1, &Qi.matrix, G2);

        // orthonormalize (Gram-Schmidt equivalent)
        kjg_gsl_matrix_QR(G2);
        kjg_gsl_matrix_QR(&Qi.matrix);

        Gswap = G2;
        G2 = G1;
        G1 = Gswap;
    }

    gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, I * L, X->m, L);
    kjg_fpca_XA(X, M, G1, &Qi.matrix);
    kjg_gsl_matrix_QR(&Qi.matrix);
    kjg_gsl_matrix_QR(Q);

    gsl_matrix_free(G1);
    gsl_matrix_free(G2);

    // PART B - compute B matrix, take SVD and return
    gsl_matrix* B = gsl_matrix_alloc(X->n, (I + 1) * L);
    kjg_fpca_XTB(X, M, Q, B);

    gsl_matrix* Vtilda = gsl_matrix_alloc(X->n, (I + 1) * L);
    gsl_matrix* Utilda = gsl_matrix_alloc((I + 1) * L, (I + 1) * L);
    gsl_vector* Stilda = gsl_vector_alloc((I + 1) * L);

    kjg_gsl_matrix_SVD(B, Utilda, Vtilda, Stilda);

    gsl_matrix_view Uk = gsl_matrix_submatrix(Utilda, 0, 0, (I + 1) * L,
            S->size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, &Uk.matrix, 0, U);

    gsl_matrix_view Vk = gsl_matrix_submatrix(Vtilda, 0, 0, X->n, S->size);
    gsl_matrix_memcpy(V, &Vk.matrix);

    gsl_vector_view Sk = gsl_vector_subvector(Stilda, 0, S->size);
    gsl_vector_memcpy(S, &Sk.vector);

    gsl_matrix_free(Q);
    gsl_matrix_free(B);
    gsl_matrix_free(Utilda);
    gsl_matrix_free(Vtilda);
    gsl_vector_free(Stilda);
}

void kjg_fpca_XTXA (
        const kjg_geno *X,
        const double* M,
        const gsl_matrix *A1,
        gsl_matrix *B,
        gsl_matrix *A2) {
    size_t i, r;                                                // row index
    uint8_t *x = malloc(sizeof(uint8_t) * X->n);                // genotypes
    double *Y = malloc(sizeof(double) * X->n * KJG_FPCA_ROWS); // normalized

    gsl_matrix_view Hmat, Xmat;

    gsl_matrix_set_zero(B);
    gsl_matrix_set_zero(A2);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        Hmat = gsl_matrix_submatrix(B, i, 0, r, B->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xmat.matrix, A1, 0,
                &Hmat.matrix);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xmat.matrix, &Hmat.matrix,
                1, A2);
    }

    free(x);
    free(Y);
}

void kjg_fpca_XA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A,
        gsl_matrix *B) {
    size_t i, r;
    uint8_t *x = malloc(sizeof(uint8_t) * X->n);
    double *Y = malloc(sizeof(double) * X->n * KJG_FPCA_ROWS);

    gsl_matrix_view Hmat, Xmat;

    gsl_matrix_set_zero(B);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        Hmat = gsl_matrix_submatrix(B, i, 0, r, B->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xmat.matrix, A, 0,
                &Hmat.matrix);
    }

    free(x);
    free(Y);
}

void kjg_fpca_XTB (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *B,
        gsl_matrix *A) {
    size_t i, r;
    uint8_t *x = malloc(sizeof(uint8_t) * X->n);
    double *Y = malloc(sizeof(double) * X->n * KJG_FPCA_ROWS);
    gsl_matrix_view Xmat;

    gsl_matrix_set_zero(A);

    for (i = 0; i < X->m; i += KJG_FPCA_ROWS) {
        r = kjg_geno_get_normalized_rows(x, Y, X, M, i, KJG_FPCA_ROWS);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        gsl_matrix_const_view Hmat = gsl_matrix_const_submatrix(B, i, 0, r,
                B->size2);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xmat.matrix, &Hmat.matrix,
                1, A);
    }

    free(x);
    free(Y);
}
