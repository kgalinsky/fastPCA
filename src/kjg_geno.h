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

// Constructor/Destructor

/** Allocates geno struct
 * @param m rows (snps)
 * @param n columns (individuals)
 */

kjg_geno* kjg_geno_alloc (size_t m, size_t n);

/** Frees geno struct
 * @param *g geno struct to free
 */

void kjg_geno_free (kjg_geno* g);

// Getter/Setter

/** Gets a row in the geno object
 * @param *g geno object
 * @param i row index
 * @param *x unpacked row
 */

void kjg_geno_get_row (const kjg_geno* g, const size_t i, uint8_t* x);

/** Sets a row in the geno object
 * @param *g geno object
 * @param i row index
 * @param *x unpacked row
 */

void kjg_geno_set_row (kjg_geno* g, const size_t i, const uint8_t* x);

// Basic functions

/** Computes mean of a row
 * @param *g geno object
 * @param i row index
 */

double kjg_geno_row_mean (const kjg_geno* g, const size_t i);

/** Compute the mean genotype of all SNPs in the geno object
 * @param *g geno object
 * @param *M array to store means
 */

void kjg_geno_row_means (const kjg_geno* g, double* M);

/** Computes the normalization lookup array.
 * @param m genotype mean
 * @param s[4] array to store the scale
 * @return success (0) or zero genotype variance error (1)
 */

int kjg_geno_normalization_lookup (const double m, double s[4]);

/** Gets a normalized row
 * @param *g geno struc
 * @param *M mean array
 * @param i row index
 * @param *y normalized genotype row
 */

void kjg_geno_get_normalized_row (
        const kjg_geno* g,
        const double* M,
        const size_t i,
        double* y);

/** Gets multiple normalized rows
 * @param *g geno struct
 * @param *M mean array
 * @param i row index
 * @param r number of rows to get
 * @param *Y normalized genotype rows
 * @return number of rows retrieved
 */

size_t kjg_geno_get_normalized_rows (
        const kjg_geno* g,
        const double* M,
        const size_t i,
        const size_t r,
        double* Y);

#endif /* KJG_GENO_H_ */
