/*
 * kjg_geno.c
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 *
 * Functions for working with genotype data.
 */

#include "kjg_geno.h"

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>

double kjg_geno_mean(const uint8_t *x, const size_t n) {
    size_t i;           // iterator
    size_t m = 0;       // good genotypes
    unsigned int t = 0; // sum of genotypes

    for (i=0; i<n; i++) {
        if (x[i] < 3) { // genotype not failure
            t += x[i];
            m++;
        }
    }

    if (m == 0) {
        return(-1);
    }

    return(((double) t) / m);
}

int kjg_geno_normalization_lookup(const double m, double s[4]) {
    double p = m / 2;
    double q = 1 - p;
    size_t i;

    s[3] = 0;

    if (p == 0 || q == 0) { // check for division by zero

        for (i=0; i<3; i++) {
            s[i] = 0; // zero out scaling array
        }

        return(1); // error
    }

    double d = sqrt(2 * p * q); // variance for Binomial(2, p) = 2pq
    for (i = 0; i < 3; i++) {
        s[i] = (i - m) / d;
    }

    return (0);
}

int kjg_geno_normalize(const uint8_t *x, double *y, const size_t n) {
    double m = kjg_geno_mean(x, n);
    int r = kjg_geno_normalize_m(m, x, y, n);
    return(r);
}

int kjg_geno_normalize_m(const double m, const uint8_t* x, double* y,
        const size_t n) {
    double s[4];
    int r = kjg_geno_normalization_lookup(m, s);
    kjg_geno_remap(s, x, y, n);
    return r;
}

void kjg_geno_remap(const double s[4], const uint8_t* x, double* y,
        const size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        y[i] = s[x[i]];
    }
}

int kjg_geno_snp_correlation(const uint8_t* x, double* C, const size_t n) {
    double m = kjg_geno_mean(x, n);

    double s[4];
    int r = kjg_geno_normalization_lookup(m, s);

    double S[4][4];
    kjg_geno_correlation_lookup(s, S);

    size_t i, j, ni;
    double l;
    uint8_t xi;
    for (i = 0; i < n; i++) {
        ni = n * i;
        xi = x[i];
        C[ni+i] += S[xi][xi];
        for (j = 0; j < i; j++) {
            l = S[xi][x[j]];
            C[ni+j] += l;
        }
    }
    return (r);
}

void kjg_geno_correlation_lookup(const double s[4], double S[4][4]) {
    size_t i, j;
    for (i = 0; i < 3; i++) {
        S[i][3] = 0;
        S[3][i] = 0;
        for (j = 0; j < 3; j++) {
            S[i][j] = S[j][i] = s[i] * s[j];
        }
    }
}

double* kjg_geno_correlation_matrix_init(const size_t n) {
    size_t i, n2 = n*n;
    double* C = malloc(sizeof(double) * n2);
    for (i = 0; i < n2; i++) {
        C[i] = 0;
    }
    return (C);
}

void kjg_geno_correlation_matrix_finish(double* C, const size_t n) {
    size_t i, j, ni;
    for (i = 1; i < n; i++) {
        ni = n*i;
        for (j = i+1; j < n; j++) {
            C[ni+j] = C[j*n+i];
        }
    }
}

void kjg_geno_pack(const uint8_t* x, uint8_t* y, size_t n) {
    size_t i, j = 0;
    uint8_t* p;
    for (i = 0; i < n; i ++) {
        if (i % 4 == 0) {
            p = y + j++;
            *p = 0;
        } else {
            *p <<= 2;
        }
        *p |= x[i];
    }
}

void kjg_geno_unpack(uint8_t* x, const uint8_t* y, size_t n) {
    size_t i, j = n/4;
    uint8_t p;

    if (n % 4 != 0) {
        p = y[j] << 2;
    }

    for (i = n; i > 0; ) {
        if (i % 4 == 0) {
            p = y[--j];
        } else {
            p >>= 2;
        }
        x[--i] = p & 3;
    }
}

void kjg_geno_free(kjg_geno* g) {
    free(g->data);
    free(g);
}

kjg_geno* kjg_geno_alloc(size_t m, size_t n) {
    kjg_geno* g = malloc(sizeof(kjg_geno));
    g->m = m;
    g->n = n;
    g->tda = (n + 3) / 4;
    g->data = malloc(sizeof(uint8_t) * g->m * g->tda);
    return (g);
}

void kjg_geno_get_row(uint8_t* x, const kjg_geno* g, const size_t i) {
    kjg_geno_unpack(x, g->data + g->tda * i, g->n);
}

void kjg_geno_set_row(const uint8_t* x, kjg_geno* g, const size_t i) {
    kjg_geno_pack(x, g->data + (g->tda * i), g->n);
}
void kjg_geno_row_means(const kjg_geno* g, double* M, uint8_t* x) {
    size_t i;
    for (i = 0; i < g->m; i++) {
        kjg_geno_get_row(x, g, i);
        M[i] = kjg_geno_mean(x, g->n);
    }
}
