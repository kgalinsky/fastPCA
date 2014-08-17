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
#include <string.h>

#include "kjg_2bit.h"

// Constructor/Destructor

kjg_geno* kjg_geno_alloc (size_t m, size_t n) {
    kjg_geno pre = { m, n, kjg_2bit_packed_tda(n) };
    kjg_geno* g = malloc(sizeof(kjg_geno));

    memcpy(g, &pre, sizeof(kjg_geno));
    g->data = malloc(sizeof(uint8_t) * g->m * g->tda);

    return (g);
}

void kjg_geno_free (kjg_geno* g) {
    free(g->data);
    free(g);
}

// Getter/Setter

void kjg_geno_get_row (const kjg_geno* g, const size_t i, uint8_t* x) {
    kjg_2bit_unpack(g->n, g->data + g->tda * i, x);
}

void kjg_geno_set_row (kjg_geno* g, const size_t i, const uint8_t* x) {
    kjg_2bit_pack(g->n, x, g->data + g->tda * i);
}

// Mean

#define MEAN(n) (((n) & 3) % 3) + ((((n) >> 2) & 3) % 3) + ((((n) >> 4) & 3) % 3) + ((((n) >> 6) & 3) % 3)

// macros for unpacking array generation
#define M1(n) MEAN(n),      MEAN(n | 1),      MEAN(n | 2),      MEAN(n | 3)
#define M2(n)   M1(n), M1(n | (1 << 2)), M1(n | (2 << 2)), M1(n | (3 << 2))
#define M3(n)   M2(n), M2(n | (1 << 4)), M2(n | (2 << 4)), M2(n | (3 << 4))

// mean lookup array
static const uint8_t MEAN_LOOKUP[256] =
    { M3(0), M3((1 << 6)), M3((2 << 6)), M3((3 << 6)) };

double kjg_geno_row_mean (const kjg_geno* g, const size_t i) {
    size_t tda = g->tda;
    uint8_t* p = g->data + (g->tda * i);
    size_t sum = 0;
    while (tda--) sum += MEAN_LOOKUP[*(p++)];
    return( ((double) sum) / g->n);
}

void kjg_geno_row_means (const kjg_geno* g, double* M) {
    size_t i;
    for (i = 0; i < g->m; i++) M[i] = kjg_geno_row_mean(g, i);
}

// Normalization

int kjg_geno_normalization_lookup (const double m, double s[4]) {
    size_t i;
    s[3] = 0;

    if (m == 0 || m == 2) {             // check for homogeneous population
        for (i = 0; i < 3; i++)
            s[i] = 0;                   // zero out scaling array
        return (1);                     // error
    }

    double p = m / 2;                   // G ~ Binomial(2, p)
    double d = sqrt(2 * p * (1 - p));   // Var(G) = 2pq

    for (i = 0; i < 3; i++)
        s[i] = (i - m) / d;

    return (0);
}

void kjg_geno_get_normalized_row (
        const kjg_geno* g,
        const double* M,
        const size_t i,
        double* y) {
    double s[4];
    kjg_geno_normalization_lookup(M[i], s);

    size_t tda = g->tda;
    uint8_t* p = g->data + (g->tda * i);

    size_t j;

    while (--tda) {
        const uint8_t* u = KJG_2BIT_UNPACK_LOOKUP[*(p++)];
        for (j = 0; j < 4; j++) *(y++) = s[*(u++)];
    }

    const uint8_t* u = KJG_2BIT_UNPACK_LOOKUP[*p];
    for (j = (g->tda - 1 ) * 4; j < i; j++) *(y++) = s[*(u++)];
}

size_t kjg_geno_get_normalized_rows (
        const kjg_geno* g,
        const double* M,
        const size_t i,
        const size_t r,
        double* Y) {
    size_t j;
    for (j = i; j < i + r && j < g->m; j++) {
        kjg_geno_get_normalized_row(g, M, j, Y);
        Y += g->n;
    }
    return (j - i);
}
