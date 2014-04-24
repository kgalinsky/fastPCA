/*
 * kjg_geno.h
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 */

#ifndef KJG_GENO_H_
#define KJG_GENO_H_

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t m;
    size_t n;
    size_t tda;
    uint8_t *data;
} kjg_geno;

/**
 * Compute the mean of the genotypes.
 *
 * @param *x array of genotypes
 * @param n number of subjects
 * @return mean of the genotypes whose assays didn't fail
 */
double kjg_geno_mean(const uint8_t *x, const size_t n);

/**
 * Normalize genotypes.
 *
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */
int kjg_geno_normalize(const uint8_t *x, double *y, const size_t n);

/**
 * Normalize genotypes when you already have the mean.
 *
 * @param m mean of the genotypes
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */
int kjg_geno_normalize_m(const double m, const uint8_t* x, double* y,
        const size_t n);

/**
 * Remap genotypes values to doubles.
 *
 * @param s array of genotype mappings
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success (0) or zero-good genos (1)
 */
void kjg_geno_remap(const double s[4], const uint8_t* x, double* y,
        const size_t n);

/**
 * Compute the normalization lookup array.
 *
 * @param m genotype mean
 * @param s[4] array to store the scale
 * @return success (0) or zero-good genos (1)
 */
int kjg_geno_normalization_lookup(const double m, double s[4]);

/**
 * Initialize GRM.
 *
 * @param n number of subjects
 * @return array
 */
double* kjg_geno_GRM_init(const size_t n);

/**
 * Add correlations for a SNP to the GRM.
 *
 * @params *x array of genotypes
 * @params *GRM genetic relationship matrix
 * @params n number of subjects
 * @return success (0) or zero-good genos (1)
 */
int kjg_geno_GRM_update(const uint8_t* x, double* GRM, const size_t n);

/**
 * Compute the GRM lookup array.
 *
 * @param s[4] normalization lookup
 * @param S[4][4] 4x4 array to store products of scaled genotypes
 * @return success (0) or zero-good genos (1)
 */
void kjg_geno_GRM_lookup(const double s[4], double S[4][4]);

/**
 * Pack genotype array
 *
 * @param *x raw genotypes
 * @param *y where to pack them
 * @param n number of genotypes
 */
void kjg_geno_pack(const uint8_t* x, uint8_t* y, size_t n);

/**
 * Unack genotype array
 *
 * @param *x where to unpack genotypes
 * @param *y packed genotypes
 * @param n number of genotypes
 */

void kjg_geno_unpack(uint8_t* x, const uint8_t* y, size_t n);

/**
 * Allocate geno struct
 *
 * @param m rows (snps)
 * @param n columns (individuals)
 */

kjg_geno* kjg_geno_alloc(size_t m, size_t n);

/**
 * Free geno struct
 *
 * @param *g geno struct to free
 */

void kjg_geno_free(kjg_geno* g);

/**
 * Set a row in the geno object
 *
 * @param *x unpacked genotype row
 * @param *g geno object
 * @param i row index
 */

void kjg_geno_set_row(const uint8_t* x, kjg_geno* g, const size_t i);

/**
 * Get a row in the geno object
 *
 * @param *x unpacked genotype row
 * @param *g geno object
 * @param i row index
 */

void kjg_geno_get_row(uint8_t* x, const kjg_geno* g, const size_t i);

/**
 * Compute the mean genotype of all SNPs in the geno object
 *
 * @param *g geno object
 * @param *M array to store means
 */

void kjg_geno_row_means(const kjg_geno* g, double* M);

#endif /* KJG_GENO_H_ */
