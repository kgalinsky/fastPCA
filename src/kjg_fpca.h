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
 * @param *M SNP means
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param L width of projection matrix
 * @param I iterations to do exponentiation
 */

void kjg_fpca (
        const kjg_geno* X,
        const double* M,
        gsl_vector* eval,
        gsl_matrix* evec,
        size_t L,
        size_t I);

#endif /* KJG_FPCA_H_ */
