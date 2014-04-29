/*
 * kjg_fpca.h
 *
 *  Created on: Apr 28, 2014
 *      Author: Kevin
 */

#ifndef KJG_FPCA_H_
#define KJG_FPCA_H_

#include <gsl/gsl_matrix.h>
#include "kjg_geno.h"

extern size_t KJG_FPCA_ROWS; // number of rows to process at once

/**
 * FastPCA blanczos step
 *
 * @param *X compressed genotype matrix (MxN)
 * @param *M array of SNP means
 * @param *G random norm matrix (NxL)
 * @param *H matrix to store product (MxIL)
 */

void kjg_fpca_blanczos(const kjg_geno* X, const double *M, gsl_matrix* G,
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

void kjg_fpca_XTXG(const kjg_geno *X, const double *M, const gsl_matrix *G1,
        gsl_matrix *H, gsl_matrix *G2);

/**
 * Multiply H = X*G
 *
 * @param *X compressed genotype matrix
 * @param *M array of SNP means
 * @param *G some matrix
 * @param *H another matrix
 */

void kjg_fpca_XG(const kjg_geno *X, const double *M, const gsl_matrix *G,
        gsl_matrix *H);

/**
 * Multiply G = XT*H
 *
 * @param *X compressed genotype matrix
 * @param *M array of SNP means
 * @param *H some matrix
 * @param *G another matrix
 */

void kjg_fpca_XTH(const kjg_geno *X, const double *M, const gsl_matrix *H,
        gsl_matrix *G);

#endif /* KJG_FPCA_H_ */
