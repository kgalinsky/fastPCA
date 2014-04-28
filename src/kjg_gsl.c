/*
 * kjg_gsl.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include "kjg_gsl.h"

void kjg_gsl_matrix_fprintf(FILE* stream, gsl_matrix* m, const char* template) {
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

void kjg_gsl_evec_fprintf(FILE* stream, gsl_vector* eval, gsl_matrix* evec,
        const char* template) {
    size_t i, j;
    fprintf(stream, "#");
    fprintf(stream, template, gsl_vector_get(eval, 0));
    for (i = 1; i < eval->size; i++) {
        fprintf(stream, "\t");
        fprintf(stream, template, gsl_vector_get(eval, i));
    }
    fprintf(stream, "\n");
    kjg_gsl_matrix_fprintf(stream, evec, template);
}

gsl_rng *kjg_gsl_rng_init() {
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    fprintf(stderr, "generator type: %s\n", gsl_rng_name(r));
    fprintf(stderr, "seed = %lu\n", gsl_rng_default_seed);

    return(r);
}

int kjg_gsl_matrix_frobenius_normalize(gsl_matrix* m) {
    double s = kjg_gsl_matrix_frobenius_norm(m);
    double d = m->size1 * m->size2;
    return (gsl_matrix_scale(m, d / s));
}

float kjg_gsl_matrix_frobenius_norm(const gsl_matrix* m) {
    if (m->size2 == m->tda) {
        gsl_vector_const_view V = gsl_vector_const_view_array(m->data, m->size1*m->size2);
        return(gsl_blas_dnrm2(&V.vector));
    }

    size_t i;
    double norm = 0;
    for (i = 0; i < m->size1; i++) {
        gsl_vector_const_view V = gsl_matrix_const_row (m, i);
        double n = gsl_blas_dnrm2(&V.vector);
        norm += n*n;
    }
    return(sqrt(norm));
}

void kjg_gsl_matrix_set_ran_ugaussian(gsl_matrix* m, const gsl_rng* r) {
    size_t i, j;
    double x, y, r2;
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j+=2) {
            do {
                /* choose x,y in uniform square (-1,-1) to (+1,+1) */
                x = -1 + 2 * gsl_rng_uniform_pos (r);
                y = -1 + 2 * gsl_rng_uniform_pos (r);

                /* see if it is in the unit circle */
                r2 = x * x + y * y;
            } while (r2 > 1.0 || r2 == 0);
            r2 = sqrt (-2.0 * log (r2) / r2);

            gsl_matrix_set(m, i, j, x*r2);
            gsl_matrix_set(m, i, j+1, y*r2);
        }
    }
}
