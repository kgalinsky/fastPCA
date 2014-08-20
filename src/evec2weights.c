#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "kjg_geno.h"
#include "kjg_genoIO.h"

#include "kjg_gsl.h"
#include "kjg_fpca.h"

void main (int argc, char **argv) {
    if (argc != 4) {
        exit(1);
    }

    FILE *fh_evec = fopen(argv[1], "r");
    kjg_genoIO *gp = kjg_genoIO_fopen(argv[2], "r");

    size_t k = 1;
    {
        char c;

        while (1) {
            c = getc(fh_evec);
            if (c == '\t') k++;
            else if (c == ' ') k++;
            else if (c == '\n') break;
        }

        fseek(fh_evec, 0, SEEK_SET);
    }

    gsl_vector *eval = gsl_vector_alloc(k);
    gsl_matrix *evec = gsl_matrix_alloc(gp->n, k);
    gsl_matrix *weights = gsl_matrix_alloc(gp->m, k);

    printf("Reading evec (%dx%d)\n", gp->n, k);
    int r = kjg_gsl_evec_fscanf(fh_evec, eval, evec);
    if (r) {
        printf("Error reading evec file\n");
        exit(1);
    }
    fclose(fh_evec);

    printf("Reading geno (%dx%d)\n", gp->m, gp->n);
    kjg_geno* X = kjg_genoIO_fread_geno(gp);
    kjg_genoIO_fclose(gp);

    double *M = malloc(sizeof(double) * X->m);

    printf("Calculating SNP allele frequencies");
    kjg_geno_row_means(X, M);

    printf("Computing weights (%dx%d)\n", X->m, k);
    kjg_fpca_XA(X, M, evec, weights);
    gsl_matrix_free(evec);
    free(M);

    printf("Scaling weights\n");
    {
        size_t i;
        for (i = 0; i < k; i++) {
            double scale = 1.0 / sqrt(X->m * gsl_vector_get(eval, i));
            gsl_vector_view V = gsl_matrix_column(weights, i);
            gsl_vector_scale(&V.vector, scale);
        }
    }

    kjg_geno_free(X);

    printf("Printing weights\n");
    FILE *fh_weight = fopen(argv[3], "w");
    kjg_gsl_evec_fprintf(fh_weight, eval, weights, "%g");

    gsl_vector_free(eval);
    gsl_matrix_free(weights);
}
