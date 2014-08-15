/** @file kjg_fpca.h
 * @brief Runs fastPCA.
 * This module also has methods to multiply a genotype matrix against the GSL
 * matrices.
 */

#ifndef KJG_FPCA_H_
#define KJG_FPCA_H_

#include <gsl/gsl_matrix.h>
#include "kjg_geno.h"

extern size_t KJG_FPCA_ROWS; // number of rows to process at once

/** Performs a fast PCA
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

/** Performs a fast SVD
 * @param *X compressed genotype matrix (MxN)
 * @param *M SNP means
 * @param *S eigenvalues
 * @param *U eigenvectors
 * @param *V eigenvectors
 * @param L width of projection matrix
 * @param I iterations to do exponentiation
 */

void kjg_SVD (
        const kjg_geno* X,
        const double* M,
        gsl_vector* S,
        gsl_matrix* U,
        gsl_matrix* V,
        size_t L,
        size_t I);

/** Multiplies B=X*A1 and A2 = XT*B = XT*X*A1
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A1 some matrix
 * @param *B intermediate matrix
 * @param *A2 next matrix
 */

void kjg_fpca_XTXA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A1,
        gsl_matrix *B,
        gsl_matrix *A2);

/** Multiplies B = X*A
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A some matrix
 * @param *B another matrix
 */

void kjg_fpca_XA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A,
        gsl_matrix *B);

/** Multiplies A = XT*B
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *B some matrix
 * @param *A another matrix
 */

void kjg_fpca_XTB (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *B,
        gsl_matrix *A);

#endif /* KJG_FPCA_H_ */
