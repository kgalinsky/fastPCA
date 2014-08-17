/*
 * kjg_GRM.h
 *
 *  Created on: Apr 24, 2014
 *      Author: Kevin
 */

#ifndef KJG_GRM_H_
#define KJG_GRM_H_

#include <stdint.h>

#include "kjg_geno.h"

typedef struct {
    const size_t n;
    double* data;
} kjg_GRM;

/** Initializes GRM
 * @param n number of subjects
 * @return array
 */
kjg_GRM* kjg_GRM_alloc (const kjg_geno* g);

/** Frees the GRM
 * @param *GRM
 */
void kjg_GRM_free (kjg_GRM* GRM);

/** Calculates the GRM
 * @params *GRM
 * @params *g genotype matrix
 * @params M genotype means
 */
void kjg_GRM_calc (kjg_GRM* GRM, const kjg_geno* g, const double* M);

/** Compute the GRM lookup array.
 * @param s[4] normalization lookup
 * @param S[4][4] 4x4 array to store products of scaled genotypes
 * @return success (0) or zero-good genos (1)
 */
void kjg_GRM_lookup (const double s[4], double S[4][4]);

#endif /* KJG_GRM_H_ */
