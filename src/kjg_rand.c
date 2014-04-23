/*
 * kjg_rand.c
 *
 *  Created on: Apr 23, 2014
 *      Author: Kevin
 */

#include <math.h>
#include <stdlib.h>

#include "kjg_rand.h"

double kjg_runif() { return((double)rand() / RAND_MAX); }

void kjg_rnorms(size_t n, double *Z) {
    size_t i;
    double x, y, r2;
    // Polar Box-Muller algorithm generates 2 at a time
    for (i = 0; i < n - 1; i += 2) {
        do {
            x = -1 + 2*kjg_runif();
            y = -1 * 2*kjg_runif();
            r2 = x*x + y*y;
        } while (r2 > 1.0 || r2 == 0);
        r2 = sqrt(-2*log(r2)/r2);
        Z[i]   = x * r2;
        Z[i+1] = y * r2;
    }

    // If n is odd
    if (i != n) {
        do {
            x = -1 + 2*kjg_runif();
            y = -1 * 2*kjg_runif();
            r2 = x*x + y*y;
        } while (r2 > 1.0 || r2 == 0);
        r2 = sqrt(-2*log(r2)/r2);
        Z[i]   = x * r2;
    }
}
