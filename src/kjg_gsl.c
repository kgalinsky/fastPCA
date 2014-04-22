/*
 * kjg_gsl.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kjg_gsl.h"

void kjg_matrix_fprintf(FILE* stream, gsl_matrix* m, const char* template) {
    size_t i, j;
    for (i = 0; i < m->size1; i++) {
        fprintf(stream, template, gsl_matrix_get(m, i, 0));
        for (j = 1; j < m->size2; j++) {
            fprintf(stream, "\t");
            fprintf(stream, template, gsl_matrix_get(m, i, j));
        }
        fprintf(stream, "\n");
    }
}

gsl_rng *kjg_rng_init() {
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    fprintf(stderr, "generator type: %s\n", gsl_rng_name(r));
    fprintf(stderr, "seed = %lu\n", gsl_rng_default_seed);

    return(r);
}

int kjg_frobenius_normalize(gsl_matrix* m) {
    double s = kjg_frobenius_norm(m);
    double d = m->size1 * m->size2;
    return (gsl_matrix_scale(m, d / s));
}

float kjg_frobenius_norm(const gsl_matrix* m) {
    size_t i, j;
    double sumsq = 0;
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            double mij = gsl_matrix_get(m, i, j);
            sumsq += mij * mij;
        }
    }
    return (sqrt(sumsq));
}

void kjg_matrix_set_ran_ugaussian(gsl_matrix* m, const gsl_rng* r) {
    size_t i, j;
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            gsl_matrix_set(m, i, j, gsl_ran_ugaussian(r));
        }
    }
}

void kjg_blanczos(
        const kjg_geno* X, const double* M, gsl_matrix* G,
        gsl_matrix* H) {
    size_t i;

    gsl_matrix* G2 = gsl_matrix_alloc(G->size1, G->size2);
    gsl_matrix* Gswap;
    gsl_matrix_view Hsub;

    for (i = 0; i < H->size2; i += G->size2) {
        Hsub = gsl_matrix_submatrix(H, 0, i, H->size1, G->size2);
        kjg_XTXG(X, M, G, &Hsub.matrix, G2);
        kjg_frobenius_normalize(&Hsub.matrix);

        Gswap = G2;
        G2 = G;
        G = Gswap;
    }

    gsl_matrix_free(G2);
}

void kjg_XTXG(
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

void kjg_XG(const kjg_geno *X, const double *M, const gsl_matrix *G,
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

void kjg_XTH(const kjg_geno *X, const double *M, const gsl_matrix *H,
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

void kjg_evec_fprintf(FILE* stream, gsl_vector* eval, gsl_matrix* evec,
        const char* template) {
    size_t i, j;
    fprintf(stream, "#");
    fprintf(stream, template, gsl_vector_get(eval, 0));
    for (i = 1; i < eval->size; i++) {
        fprintf(stream, "\t");
        fprintf(stream, template, gsl_vector_get(eval, i));
    }
    fprintf(stream, "\n");
    kjg_matrix_fprintf(stream, evec, template);
}
