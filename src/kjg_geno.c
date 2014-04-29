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
    size_t i;           // index
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
    size_t i;           // index

    s[3] = 0;

    if (m == 0 || m == 2) { // check for homogeneous population
        for (i=0; i<3; i++) {
            s[i] = 0; // zero out scaling array
        }
        return(1); // error
    }

    double p = m/2;         // G ~ Binomial(2, p)
    double q = 1-p;
    double d = sqrt(2*p*q); // Var(G) = 2pq

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

void kjg_geno_get_normalized_row(uint8_t* x, double* y,
        const kjg_geno* g, const double* M, const size_t i) {
    kjg_geno_get_row(x, g, i);
    kjg_geno_normalize_m(M[i], x, y, g->n);
}

size_t kjg_geno_get_normalized_rows(uint8_t* x, double* Y,
        const kjg_geno* g, const double* M, const size_t i, const size_t r) {
    size_t j;
    double *y = Y;
    for (j = i; j < i + r && j < g->m; j++) {
        kjg_geno_get_normalized_row(x, y, g, M, j);
        y += g->n;
    }
    return(j-i);
}

void kjg_geno_row_means(const kjg_geno* g, double* M) {
    size_t i;
    uint8_t *x = malloc(sizeof(uint8_t)*g->n);

    for (i = 0; i < g->m; i++) {
        kjg_geno_get_row(x, g, i);
        M[i] = kjg_geno_mean(x, g->n);
    }

    free(x);
}
