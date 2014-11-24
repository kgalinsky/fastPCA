/**
 * @file kjg_gsl.h
 * @brief Augment GSL functions
 */

#ifndef KJG_GSL_H_
#define KJG_GSL_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

/**
 * Prints the matrix tab-delimited
 * @param *stream output file pointer
 * @param *m gsl_matrix to print
 * @param *template character template for fprintf
 */

void
kjg_gsl_matrix_fprintf (FILE* stream, gsl_matrix* m, const char* template);

/**
 * Prints the eigenvalues and then eigenvectors below
 * @param *stream output file pointer
 * @param *eval eigenvalues
 * @param *evec eigenvectors
 * @param *template character template for fprintf */

void
kjg_gsl_evec_fprintf (FILE* stream, gsl_vector* eval, gsl_matrix* evec,
                      const char* template);

/**
 * Reads a matrix
 * @param *stream input file pointer
 * @param *m matrix to store
 */

void
kjg_gsl_matrix_fscanf (FILE* stream, gsl_matrix* m);

/**
 * Reads an evec
 * @param *stream input file pointer
 * @param *eval eigenvalues vector
 * @param *evec eigenvectors matrix
 */

int
kjg_gsl_evec_fscanf (FILE* stream, gsl_vector* eval, gsl_matrix* evec);

/**
 * Initializes random number generation.
 */

gsl_rng *
kjg_gsl_rng_init ();

/**
 * Initializes the matrix with random unit gaussians
 * @param *m matrix to be set
 * @param *r random number generator
 */

void
kjg_gsl_ran_ugaussian_pair (const gsl_rng* r, double x[2]);

/**
 * Fills a matrix with unit Gaussian random variates
 * @param *r random number generator
 * @param *m matrix to be filled
 */

void
kjg_gsl_ran_ugaussian_matrix (const gsl_rng* r, gsl_matrix* m);

/**
 * Performs the QR decomposition on the matrix and return Q in the matrix
 * @param *m matrix to orthogonalize
 */

void
kjg_gsl_matrix_QR (gsl_matrix* m);

#endif /* KJG_GSL_H_ */
