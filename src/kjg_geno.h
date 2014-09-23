/**
 * @file kjg_geno.h
 * @brief Data structure and methods to store genotype data
 */

#ifndef KJG_GENO_H_
#define KJG_GENO_H_

#include <stddef.h>
#include <stdint.h>

// Data structure for genotype data

typedef struct {
    const size_t m;     // number of SNPs
    const size_t n;     // number of individuals
    const size_t tda;   // width of a packed SNP row
    uint8_t *data;      // packed genotype data
    double *af;         // allele frequencies
    double *norm;       // normalization tables - store to perform random access
} kjg_geno;

// Packing/unpacking macros
#define KJG_GENO_PACK(a, b, c, d) \
    ( (((a) & 3) << 0) | \
      (((b) & 3) << 2) | \
      (((c) & 3) << 4) | \
      (((d) & 3) << 6) )

#define KJG_GENO_UNPACK_I(p, i) (((p) >> ((i) * 2)) & 3)

#define KJG_GENO_UNPACK(p) \
    { KJG_GENO_UNPACK_I(p, 0), \
      KJG_GENO_UNPACK_I(p, 1), \
      KJG_GENO_UNPACK_I(p, 2), \
      KJG_GENO_UNPACK_I(p, 3) }

#define KJG_GENO_ALT(u) ((u) > 1 ? (u) - 1 : 0)

#define KJG_GENO_SUM_ALT(p) \
    ( KJG_GENO_ALT(KJG_GENO_UNPACK_I(p, 0)) + \
      KJG_GENO_ALT(KJG_GENO_UNPACK_I(p, 1)) + \
      KJG_GENO_ALT(KJG_GENO_UNPACK_I(p, 2)) + \
      KJG_GENO_ALT(KJG_GENO_UNPACK_I(p, 3)) )

#define KJG_GENO_PRESENT(u) ((u) != 1)

#define KJG_GENO_COUNT(p) \
    ( KJG_GENO_PRESENT(KJG_GENO_UNPACK_I(p, 0)) + \
      KJG_GENO_PRESENT(KJG_GENO_UNPACK_I(p, 1)) + \
      KJG_GENO_PRESENT(KJG_GENO_UNPACK_I(p, 2)) + \
      KJG_GENO_PRESENT(KJG_GENO_UNPACK_I(p, 3)) )

extern const uint8_t KJG_GENO_PACK_LOOKUP[4][4][4][4];
extern const uint8_t KJG_GENO_UNPACK_LOOKUP[256][4];
extern const uint8_t KJG_GENO_SUM_ALT_LOOKUP[256];
extern const uint8_t KJG_GENO_COUNT_LOOKUP[256];

// Functional

void kjg_geno_pack (const size_t n, const uint8_t* u, uint8_t* p);
void kjg_geno_unpack (const size_t n, const uint8_t* p, uint8_t* u);
size_t kjg_geno_sum_alt (const size_t n, const uint8_t* p);
size_t kjg_geno_count (const size_t n, const uint8_t* p);
double kjg_geno_af (const size_t n, const uint8_t* p);

/**
 * Computes the normalization lookup array.
 * @param p alt allele frequency
 * @param s[4] array to store the scale
 * @return success (0) or zero genotype variance error (1)
 */

int kjg_geno_norm (const double p, double s[4]);

// Constructor/Destructor

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

// Getter/Setter

/**
 * Gets a row in the geno object
 * @param *g geno object
 * @param i row index
 * @param *x unpacked row
 */

void kjg_geno_get_row (const kjg_geno* g, const size_t i, uint8_t* x);

/**
 * Gets a normalized row
 * @param *g geno struct
 * @param i row index
 * @param *y normalized genotype row
 */

void kjg_geno_get_row_normalized (const kjg_geno* g, const size_t i, double* y);

size_t kjg_geno_get_rows_normalized (const kjg_geno* g, const size_t i, const size_t r, double* Y);


/**
 * Sets a row in the geno object
 * @param *g geno object
 * @param i row index
 * @param *x unpacked row
 */

void kjg_geno_set_row (kjg_geno* g, const size_t i, const uint8_t* x);

void kjg_geno_set_af (kjg_geno* g);
void kjg_geno_set_norm (kjg_geno* g);

#endif /* KJG_GENO_H_ */
