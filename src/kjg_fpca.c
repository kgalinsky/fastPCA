/*
 * kjg_fpca.c
 *
 *  Created on: Apr 28, 2014
 *      Author: Kevin
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

#include "kjg_fpca.h"
#include "kjg_geno.h"
#include "kjg_geno_gsl.h"
#include "kjg_gsl.h"
#include "glue.h"

void
kjg_fpca (const kjg_geno* X, gsl_vector* eval, gsl_matrix* evec, size_t L,
          size_t I)
{

  if (evec->size1 != X->n)
    exit (1);
  if (eval->size != evec->size2)
    exit (1);
  if (eval->size >= L)
    exit (1);
  if (I == 0)
    exit (1);

  // PART A - compute Q such that X ~ Q * (Q^T) * X
  gsl_matrix* Q = kjg_fpca_subspace_iteration_blanczos (X, L, I);
//    glue_dgesdd('O', Q, NULL, NULL, NULL);
  glue_dgesvd ('O', 'N', Q, NULL, NULL, NULL, NULL);

  // kjg_gsl_matrix_QR(Q); // QR decomposition is less accurate than SVD

  // PART B - compute B matrix, take SVD and return
  gsl_matrix* B = gsl_matrix_alloc (X->n, (I + 1) * L);
  kjg_geno_gsl_XTB (X, Q, B);

  gsl_vector* Sl = gsl_vector_alloc (B->size2);
  gsl_matrix* Vl = B;
//    glue_dgesdd('O', B, Sl, NULL, NULL);
  glue_dgesvd ('O', 'N', B, Sl, NULL, NULL, NULL);

  gsl_matrix_view Vk = gsl_matrix_submatrix (Vl, 0, 0, X->n, eval->size);
  gsl_matrix_memcpy (evec, &Vk.matrix);

  gsl_vector_view Sk = gsl_vector_subvector (Sl, 0, eval->size);
  gsl_vector_mul (&Sk.vector, &Sk.vector);
  gsl_vector_scale (&Sk.vector, 1.0 / X->m);
  gsl_vector_memcpy (eval, &Sk.vector);

  gsl_matrix_free (Q);
  gsl_matrix_free (B);
  gsl_vector_free (Sl);
}

gsl_matrix*
kjg_fpca_subspace_iteration_blanczos (const kjg_geno* X, size_t l, size_t q)
{
  size_t m = X->m;
  size_t n = X->n;
  size_t ql = l * (q + 1);

  gsl_matrix* G1 = gsl_matrix_alloc (n, l);
  gsl_matrix* G2 = gsl_matrix_alloc (n, l);
  gsl_matrix* Gs;

  gsl_matrix* Q = gsl_matrix_alloc (m, ql);
  gsl_matrix_view Qi;

  gsl_rng *r = kjg_gsl_rng_init ();
  kjg_gsl_ran_ugaussian_matrix (r, G1);
  gsl_rng_free (r);

  kjg_gsl_matrix_QR (G1);

  size_t i;
  for (i = 0; i < q; i++)
    {
      Qi = gsl_matrix_submatrix (Q, 0, i * l, m, l);

      // do the multiplication
      kjg_geno_gsl_XTXA (X, G1, &Qi.matrix, G2);

      // orthonormalize (Gram-Schmidt equivalent)
      kjg_gsl_matrix_QR (G2);

      Gs = G2;
      G2 = G1;
      G1 = Gs;
    }

  Qi = gsl_matrix_submatrix (Q, 0, q * l, X->m, l);
  kjg_geno_gsl_XA (X, G1, &Qi.matrix);

  gsl_matrix_free (G1);
  gsl_matrix_free (G2);

  return (Q);
}
