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

    // Box-Muller algorithm generates 2 at a time
    for (i = 0; i < n - 1; i += 2) {
        double R = sqrt(-2 * log(kjg_runif()));
        double T = 2 * M_PI * kjg_runif();
        Z[i]   = R*cos(T);
        Z[i+1] = R*sin(T);
    }

    // If n is odd
    if (i != n) {
        double R = sqrt(-2 * log(kjg_runif()));
        double T = 2 * M_PI * kjg_runif();
        Z[i]   = R*cos(T);
    }
}
