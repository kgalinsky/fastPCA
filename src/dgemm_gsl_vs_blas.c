/*
 * foo.c
 *
 *  Created on: Apr 28, 2014
 *      Author: Kevin
 */

#include <stdio.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
//#include <cblas.h>
#include "kjg_gsl.h"

void print_time () {
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    printf("CPU: %d.%09d\t", t.tv_sec, t.tv_nsec);
    clock_gettime(CLOCK_REALTIME, &t);
    printf("Real: %d.%09d\n", t.tv_sec, t.tv_nsec);
}

void main () {
    size_t m=2000, n=2000, k=1000;

    gsl_rng *r = kjg_gsl_rng_init();

    gsl_matrix *X = gsl_matrix_alloc(m, k);
    gsl_matrix *Y = gsl_matrix_alloc(k, n);
    gsl_matrix *A = gsl_matrix_alloc(m, n);
    gsl_matrix *B = gsl_matrix_alloc(m, n);

    kjg_gsl_matrix_set_ran_ugaussian(X, r);
    kjg_gsl_matrix_set_ran_ugaussian(Y, r);

    print_time();
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, X, Y, 0, A);
    print_time();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, 1.0, X->data, k, Y->data, n, 0.0, B->data, n);
    print_time();
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, X, Y, 0, A);
    print_time();
}
