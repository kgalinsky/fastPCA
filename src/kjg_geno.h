/**
 * @file kjg_geno.h
 * @brief Data structure and methods to store genotype data
 */

#ifndef KJG_GENO_H_
#define KJG_GENO_H_

#include <stddef.h>
#include <stdint.h>

/** Data structure for genotype data */
typedef struct {
    const size_t m;     // number of SNPs
    const size_t n;     // number of individuals
    const size_t tda;   // width of a packed SNP row
    uint8_t *data;      // packed genotype data
} kjg_geno;

/**
 * Computes the mean of the genotypes.
 * @param *x array of genotypes
 * @param n number of subjects
 * @return mean of the genotypes whose assays didn't fail
 */

double kjg_geno_mean (const uint8_t *x, const size_t n);

/**
 * Normalizes genotypes.
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */

int kjg_geno_normalize (const uint8_t *x, double *y, const size_t n);
/**
 * Normalizes genotypes when you already have the mean.
 * @param m mean of the genotypes
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success or failure (-1)
 */

int kjg_geno_normalize_m (
        const double m,
        const uint8_t* x,
        double* y,
        const size_t n);

/**
 * Remaps genotypes values to doubles.
 * @param s array of genotype mappings
 * @param *x array of genotypes
 * @param *y array to put scaled genotypes
 * @param n number of subjects
 * @return success (0) or zero-good genos (1)
 */

void kjg_geno_remap (
        const double s[4],
        const uint8_t* x,
        double* y,
        const size_t n);

/**
 * Computes the normalization lookup array.
 * @param m genotype mean
 * @param s[4] array to store the scale
 * @return success (0) or zero-good genos (1)
 */

int kjg_geno_normalization_lookup (const double m, double s[4]);

/**
 * Packs genotype array
 * @param *x raw genotypes
 * @param *y where to pack them
 * @param n number of genotypes
 */

void kjg_geno_pack (const uint8_t* x, uint8_t* y, size_t n);

/**
 * Unacks genotype array
 * @param *x where to unpack genotypes
 * @param *y packed genotypes
 * @param n number of genotypes
 */

void kjg_geno_unpack (uint8_t* x, const uint8_t* y, size_t n);

/**
 * Allocates geno struct
 * @param m rows (snps)
 * @param n columns (individuals)
 */

kjg_geno* kjg_geno_alloc (size_t m, size_t n);

/**
 * Frees geno struct
 * @param *g geno struct to free
 */

void kjg_geno_free (kjg_geno* g);

/**
 * Sets a row in the geno object
 * @param *x unpacked genotype row
 * @param *g geno object
 * @param i row index
 */

void kjg_geno_set_row (const uint8_t* x, kjg_geno* g, const size_t i);

/**
 * Gets a row in the geno object
 * @param *x unpacked genotype row
 * @param *g geno object
 * @param i row index
 */

void kjg_geno_get_row (uint8_t* x, const kjg_geno* g, const size_t i);

/**
 * Gets a row and normalize it.
 * @param *x unpacked genotype row
 * @param *y normalized genotype row
 * @param *g geno struct
 * @param *M mean array
 * @param i row index
 */

void kjg_geno_get_normalized_row (
        uint8_t* x,
        double* y,
        const kjg_geno* g,
        const double* M,
        const size_t i);

/**
 * Gets multiple normalized rows
 * @param *x unpacked genotype row (will hold last one)
 * @param *Y normalized genotype rows
 * @param *g geno struct
 * @param *M mean array
 * @param i row index
 * @param r number of rows to get
 * @return number of rows retrieved
 */

size_t kjg_geno_get_normalized_rows (
        uint8_t* x,
        double* Y,
        const kjg_geno* g,
        const double* M,
        const size_t i,
        const size_t r);

/**
 * Compute the mean genotype of all SNPs in the geno object
 * @param *g geno object
 * @param *M array to store means
 */

void kjg_geno_row_means (const kjg_geno* g, double* M);

#endif /* KJG_GENO_H_ */
