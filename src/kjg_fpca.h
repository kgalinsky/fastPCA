/**
 * @file kjg_fpca.h
 * @brief Runs fastPCA.
 * This module also has methods to multiply a genotype matrix against the GSL
 * matrices.
 */

#ifndef KJG_FPCA_H_
#define KJG_FPCA_H_

#include <gsl/gsl_matrix.h>
#include "kjg_geno.h"

extern size_t KJG_FPCA_ROWS; // number of rows to process at once

/**
 * Performs a fast PCA
 * @param *X compressed genotype matrix (MxN)
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param L width of projection matrix
 * @param I iterations to do exponentiation
 */

void
kjg_fpca (const kjg_geno* X, gsl_vector* eval, gsl_matrix* evec, size_t L,
          size_t I);

/**
 * Perform a blanczos subspace iteration combining Rokhlin 2009 and Halko 2011.
 * @param X compressed genotype matrix (MxN)
 * @param l width of projection matrix
 * @param q iterations to do exponentiation
 * @return Q matrix
 */

gsl_matrix*
kjg_fpca_subspace_iteration_blanczos (const kjg_geno* X, size_t l, size_t q);

#endif /* KJG_FPCA_H_ */
