/*
 * kjg_gsl.h
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#ifndef KJG_GSL_H_
#define KJG_GSL_H_

#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

#include "kjg_geno.h"

/**
 * Prints the matrix tab-delimited
 *
 * @param *stream file pointer to print output
 * @param *m gsl_matrix to print
 * @param *template character template for fprintf
 */

void kjg_matrix_fprintf(FILE* stream, gsl_matrix* m, const char* template);

/**
 * Print the eigenvalues and then eigenvectors below
 *
 * @param *stream file pointer to print output
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param *template character template for fprintf */

void kjg_evec_fprintf(FILE* stream, gsl_vector* eval, gsl_matrix* evec, const char* template);

/**
 * Initialize random number generation.
 */

gsl_rng *kjg_rng_init();

/**
 * Initialize the matrix with random unit gaussians
 *
 * @param *m matrix to be set
 * @param *r random number generator
 */

void kjg_matrix_set_ran_ugaussian(gsl_matrix* m, const gsl_rng* r);

/**
 * Normalize the matrix so the frobenius norm is M*N
 *
 * @param *m matrix to normalize
 * @return if error
 */

int kjg_frobenius_normalize(gsl_matrix* m);

/**
 * Calculate the frobenius norm
 *
 * @param *m matrix to find norm of
 * @return norm
 */

float kjg_frobenius_norm(const gsl_matrix* m);

/**
 * FastPCA blanczos step
 *
 * @param *X compressed genotype matrix (MxN)
 * @param *M array of SNP means
 * @param *G random norm matrix (NxL)
 * @param *H matrix to store product (MxIL)
 */

void kjg_blanczos(const kjg_geno* X, const double *M, gsl_matrix* G,
        gsl_matrix* H);

/**
 * Multiply G2 = XT*H = XT*X*G1
 *
 * @param *X compressed genotype matrix
 * @param *M array of SNP means
 * @param *G1 some matrix
 * @param *H intermediate matrix
 * @param *G2 next matrix
 */

void kjg_XTXG(const kjg_geno *X, const double *M, const gsl_matrix *G1,
        gsl_matrix *H, gsl_matrix *G2);

/**
 * Multiply H = X*G
 *
 * @param *X compressed genotype matrix
 * @param *M array of SNP means
 * @param *G some matrix
 * @param *H another matrix
 */

void kjg_XG(const kjg_geno *X, const double *M, const gsl_matrix *G,
        gsl_matrix *H);

/**
 * Multiply G = XT*H
 *
 * @param *X compressed genotype matrix
 * @param *M array of SNP means
 * @param *H some matrix
 * @param *G another matrix
 */

void kjg_XTH(const kjg_geno *X, const double *M, const gsl_matrix *H,
        gsl_matrix *G);

#endif /* KJG_GSL_H_ */
