/*
 * kjg_GRM.c
 *
 *  Created on: Apr 24, 2014
 *      Author: Kevin
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "kjg_GRM.h"
#include "kjg_geno.h"

kjg_GRM* kjg_GRM_alloc (const kjg_geno* g) {
    size_t n = g->n;
    kjg_GRM pre = { g->n };
    kjg_GRM* GRM = malloc(sizeof(kjg_GRM));

    memcpy(GRM, &pre, sizeof(kjg_GRM));
    GRM->data = calloc(n * (n + 1) / 2, sizeof(double));

    return (GRM);
}

void kjg_GRM_free (kjg_GRM* GRM) {
    free(GRM->data);
    free(GRM);
}


void kjg_GRM_calc (kjg_GRM* GRM, const kjg_geno* g, const double* M) {
    size_t i, j, k;
    double* data;
    uint8_t* x = malloc(g->n * sizeof(uint8_t));

    for (i = 0; i < g->m; i++) {
        double s[4];
        int r = kjg_geno_normalization_lookup(M[i], s);
        if (r == 1) continue;

        kjg_geno_get_row(g, i, x);

        double S[4][4];
        kjg_GRM_lookup(s, S);

        data = GRM->data;
        for (j = 0; j < GRM->n; j++) {
            for (k = j; k < GRM->n; k++) {
                *(data++) += S[x[j]][x[k]];
            }
        }
    }

    free(x);
}

void kjg_GRM_lookup (const double s[4], double S[4][4]) {
    size_t i, j;
    for (i = 0; i < 3; i++) {
        S[i][3] = 0;
        S[3][i] = 0;
        for (j = 0; j < 3; j++) {
            S[i][j] = S[j][i] = s[i] * s[j];
        }
    }
}
