/*
 * kjg_fpca.c
 *
 *  Created on: Apr 28, 2014
 *      Author: Kevin
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

#include "kjg_fpca.h"
#include "kjg_geno.h"
#include "kjg_geno_gsl.h"
#include "kjg_gsl.h"

void kjg_fpca (
        const kjg_geno* X,
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
    kjg_gsl_ran_ugaussian_matrix(r, G1);
    gsl_rng_free(r);

    size_t i;
    for (i = 0; i < I; i++) {
        gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, i * L, X->m, L);

        // do the multiplication
        kjg_geno_gsl_XTXA(X, G1, &Qi.matrix, G2);

        // orthonormalize (Gram-Schmidt equivalent)
        kjg_gsl_matrix_QR(G2);

        Gswap = G2;
        G2 = G1;
        G1 = Gswap;
    }

    gsl_matrix_view Qi = gsl_matrix_submatrix(Q, 0, I * L, X->m, L);
    kjg_geno_gsl_XA(X, G1, &Qi.matrix);

    {
        gsl_matrix* V = gsl_matrix_alloc(Q->size2, Q->size2);
        gsl_vector* S = gsl_vector_alloc(Q->size2);

        kjg_gsl_SVD(Q, V, S);

        gsl_matrix_free(V);
        gsl_vector_free(S);
    }

    // kjg_gsl_matrix_QR(Q); // QR decomposition is less accurate than SVD

    gsl_matrix_free(G1);
    gsl_matrix_free(G2);

    // PART B - compute B matrix, take SVD and return
    gsl_matrix* B = gsl_matrix_alloc(X->n, (I + 1) * L);
    kjg_geno_gsl_XTB(X, Q, B);

    gsl_matrix* Utilda = gsl_matrix_alloc((I + 1) * L, (I + 1) * L);
    gsl_vector* Stilda = gsl_vector_alloc((I + 1) * L);

    kjg_gsl_SVD(B, Utilda, Stilda);

    gsl_matrix_view Vk = gsl_matrix_submatrix(B, 0, 0, X->n, eval->size);
    gsl_matrix_memcpy(evec, &Vk.matrix);

    gsl_vector_view Sk = gsl_vector_subvector(Stilda, 0, eval->size);
    gsl_vector_mul(&Sk.vector, &Sk.vector);
    gsl_vector_scale(&Sk.vector, 1.0 / X->m);
    gsl_vector_memcpy(eval, &Sk.vector);

    gsl_matrix_free(Q);
    gsl_matrix_free(B);
    gsl_matrix_free(Utilda);
    gsl_vector_free(Stilda);
}

