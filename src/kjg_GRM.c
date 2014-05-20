/*
 * kjg_GRM.c
 *
 *  Created on: Apr 24, 2014
 *      Author: Kevin
 */

#include <stdint.h>
#include <stdlib.h>

#include "kjg_GRM.h"
#include "kjg_geno.h"

double* kjg_GRM_init(const size_t n) {
    double* GRM = calloc(n*(n+1)/2, sizeof(double));
    return(GRM);
}

int kjg_GRM_update(const uint8_t* x, double* GRM, const size_t n) {
    double m = kjg_geno_mean(x, n);

    double s[4];
    int r = kjg_geno_normalization_lookup(m, s);

    double S[4][4];
    kjg_GRM_lookup(s, S);

    size_t i, j, k=0;
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            GRM[k++] += S[x[i]][x[j]];
        }
    }
    return (r);
}

void kjg_GRM_lookup(const double s[4], double S[4][4]) {
    size_t i, j;
    for (i = 0; i < 3; i++) {
        S[i][3] = 0;
        S[3][i] = 0;
        for (j = 0; j < 3; j++) {
            S[i][j] = S[j][i] = s[i] * s[j];
        }
    }
}
