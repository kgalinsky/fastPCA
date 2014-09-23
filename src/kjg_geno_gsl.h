#ifndef KJG_GENO_GSL_H_
#define KJG_GENO_GSL_H_

/**
 * Multiplies B=X*A1 and A2 = XT*B = XT*X*A1
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A1 some matrix
 * @param *B intermediate matrix
 * @param *A2 next matrix
 */

void kjg_geno_gsl_XTXA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A1,
        gsl_matrix *B,
        gsl_matrix *A2);

/**
 * Multiplies B = X*A
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A some matrix
 * @param *B another matrix
 */

void kjg_geno_gsl_XA (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *A,
        gsl_matrix *B);

/**
 * Multiplies A = XT*B
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *B some matrix
 * @param *A another matrix
 */

void kjg_geno_gsl_XTB (
        const kjg_geno *X,
        const double *M,
        const gsl_matrix *B,
        gsl_matrix *A);

#endif /* KJG_GENO_GSL_H_ */
