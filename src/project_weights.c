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
    if (argc != 5) {
        exit(1);
    }

    FILE *fh_weight = fopen(argv[1], "r");
    kjg_genoIO* gp1 = kjg_genoIO_fopen(argv[2], "r");
    kjg_genoIO* gp2 = kjg_genoIO_fopen(argv[3], "r");

    if (gp1->m != gp2->m) exit(1);

    size_t k = 1;
    {
        char c;

        while (1) {
            c = getc(fh_weight);
            if (c == '\t') k++;
            else if (c == ' ') k++;
            else if (c == '\n') break;
        }

        fseek(fh_weight, 0, SEEK_SET);
    }

    kjg_geno *X;

    gsl_vector *eval = gsl_vector_alloc(k);
    gsl_matrix *weights = gsl_matrix_alloc(gp1->m, k);
    printf("Reading weights (%dx%d)\n", gp1->m, k);
    int r = kjg_gsl_evec_fscanf(fh_weight, eval, weights);
    if (r) {
        printf("Error reading weights file\n");
        exit(1);
    }
    fclose(fh_weight);

    printf("Reading geno1 (%dx%d)\n", gp1->m, gp1->n);
    X = kjg_genoIO_fread_geno(gp1);
    kjg_genoIO_fclose(gp1);

    double *M = malloc(sizeof(double) * X->m);
    printf("Computing means\n");
    kjg_geno_row_means(X, M);
    kjg_geno_free(X);
    free(M);

    printf("Reading geno2 (%dx%d)\n", gp2->m, gp2->n);
    X = kjg_genoIO_fread_geno(gp2);
    kjg_genoIO_fclose(gp2);

    gsl_matrix *evec = gsl_matrix_alloc(X->n, k);
    printf("Computing evec (%dx%d)\n", X->m, k);
    kjg_fpca_XTB(X, M, weights, evec);
    gsl_matrix_free(weights);

    printf("Scaling evec\n");
    {
        size_t i;
        for (i = 0; i < k; i++) {
            double scale = 1.0 / sqrt(X->m * gsl_vector_get(eval, i));
            gsl_vector_view V = gsl_matrix_column(evec, i);
            gsl_vector_scale(&V.vector, scale);
        }
    }

    kjg_geno_free(X);

    printf("Printing evec\n");
    FILE *fh_evec = fopen(argv[4], "w");
    kjg_gsl_evec_fprintf(fh_evec, eval, evec, "%g");

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}
