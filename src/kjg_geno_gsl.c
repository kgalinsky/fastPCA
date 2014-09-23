#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "kjg_geno.h"
#include "kjg_geno_gsl.h"

size_t KJG_GENO_GSL_ROWS = 256;

void kjg_geno_gsl_XTXA (
        const kjg_geno *X,
        const double* M,
        const gsl_matrix *A1,
        gsl_matrix *B,
        gsl_matrix *A2) {
    size_t i, r;                                                // row index
    double *Y = malloc(sizeof(double) * X->n * KJG_GENO_GSL_ROWS);  // normalized

    gsl_matrix_view Bi, Xi;

    gsl_matrix_set_zero(A2);

    for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS) {
        r = kjg_geno_get_normalized_rows(X, M, i, KJG_GENO_GSL_ROWS, Y);
        Xi = gsl_matrix_view_array(Y, r, X->n);
        Bi = gsl_matrix_submatrix(B, i, 0, r, B->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xi.matrix, A1, 0,
                &Bi.matrix);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xi.matrix, &Bi.matrix, 1,
                A2);
    }

    free(Y);
}

void kjg_geno_gsl_XA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A,
        gsl_matrix *B) {
    size_t i, r;
    double *Y = malloc(sizeof(double) * X->n * KJG_GENO_GSL_ROWS);

    gsl_matrix_view Hmat, Xmat;

    gsl_matrix_set_zero(B);

    for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS) {
        r = kjg_geno_get_normalized_rows(X, M, i, KJG_GENO_GSL_ROWS, Y);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        Hmat = gsl_matrix_submatrix(B, i, 0, r, B->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &Xmat.matrix, A, 0,
                &Hmat.matrix);
    }

    free(Y);
}

void kjg_geno_gsl_XTB (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *B,
        gsl_matrix *A) {
    size_t i, r;
    double *Y = malloc(sizeof(double) * X->n * KJG_GENO_GSL_ROWS);
    gsl_matrix_view Xmat;

    gsl_matrix_set_zero(A);

    for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS) {
        r = kjg_geno_get_normalized_rows(X, M, i, KJG_GENO_GSL_ROWS, Y);
        Xmat = gsl_matrix_view_array(Y, r, X->n);
        gsl_matrix_const_view Hmat = gsl_matrix_const_submatrix(B, i, 0, r,
                B->size2);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &Xmat.matrix, &Hmat.matrix,
                1, A);
    }

    free(Y);
}
