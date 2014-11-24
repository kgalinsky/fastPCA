/*
 * kjg_geno_rand.h
 *
 *  Created on: Nov 23, 2014
 *      Author: Kevin Galinsky
 */

#ifndef KJG_GENO_RAND_H_
#define KJG_GENO_RAND_H_

#include <gsl/gsl_rng.h>

#include "kjg_geno.h"

kjg_geno*
kjg_geno_rand_star (gsl_rng* r, const size_t M, const size_t N, const size_t P,
                    const double FST, const double *MAF);

double
kjg_geno_rand_anc (gsl_rng* r, const double* MAF);

void
kjg_geno_rand_star_AF (gsl_rng* r, double* AF, const double anc,
                       const double FST, const size_t P);

void
kjg_geno_rand_row (gsl_rng* r, uint8_t* x, const size_t N, const size_t P,
                   const double* AF);

void
kjg_geno_rand_ld_row (gsl_rng* r, uint8_t* x, const uint8_t *y,
                      const double cor, const size_t N, const size_t P,
                      const double* AF);

#endif /* KJG_GENO_RAND_H_ */
