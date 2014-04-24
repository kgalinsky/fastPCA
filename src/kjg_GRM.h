/*
 * kjg_GRM.h
 *
 *  Created on: Apr 24, 2014
 *      Author: Kevin
 */

#ifndef KJG_GRM_H_
#define KJG_GRM_H_

#include <stdint.h>

/**
 * Initialize GRM.
 *
 * @param n number of subjects
 * @return array
 */
double* kjg_GRM_init(const size_t n);

/**
 * Add correlations for a SNP to the GRM.
 *
 * @params *x array of genotypes
 * @params *GRM genetic relationship matrix
 * @params n number of subjects
 * @return success (0) or zero-good genos (1)
 */
int kjg_GRM_update(const uint8_t* x, double* GRM, const size_t n);

/**
 * Compute the GRM lookup array.
 *
 * @param s[4] normalization lookup
 * @param S[4][4] 4x4 array to store products of scaled genotypes
 * @return success (0) or zero-good genos (1)
 */
void kjg_GRM_lookup(const double s[4], double S[4][4]);

#endif /* KJG_GRM_H_ */
