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
    FILE *fh_geno1 = fopen(argv[2], "r");
    FILE *fh_geno2 = fopen(argv[3], "r");

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

    size_t n1 = kjg_genoIO_num_ind(fh_geno1);
    size_t n2 = kjg_genoIO_num_ind(fh_geno2);
    size_t m = kjg_genoIO_num_snp(fh_geno1, n1);

    if (m != kjg_genoIO_num_snp(fh_geno2, n2)) {
        exit(1);
    }

    gsl_vector *eval = gsl_vector_alloc(k);
    gsl_matrix *evec = gsl_matrix_alloc(n2, k);
    gsl_matrix *weights = gsl_matrix_alloc(m, k);

    kjg_geno *X;
    double *M = malloc(sizeof(double) * m);

    printf("Reading weights (%dx%d)\n", m, k);
    int r = kjg_gsl_evec_fscanf(fh_weight, eval, weights);
    if (r) {
        printf("Error reading weights file\n");
        exit(1);
    }
    fclose(fh_weight);

    printf("Reading geno1 (%dx%d)\n", m, n1);
    X = kjg_geno_alloc(m, n1);
    kjg_genoIO_fread_geno(X, fh_geno1);
    fclose(fh_geno1);

    printf("Computing means\n");
    kjg_geno_row_means(X, M);
    kjg_geno_free(X);

    printf("Reading geno2 (%dx%d)\n", m, n2);
    X = kjg_geno_alloc(m, n2);
    kjg_genoIO_fread_geno(X, fh_geno2);
    fclose(fh_geno2);

    printf("Computing evec (%dx%d)\n", m, k);
    kjg_fpca_XTH(X, M, weights, evec);

    printf("Scaling evec\n");
    {
        size_t i;
        for (i = 0; i < k; i++) {
            double scale = 1.0 / sqrt(m * gsl_vector_get(eval, i));
            gsl_vector_view V = gsl_matrix_column(evec, i);
            gsl_vector_scale(&V.vector, scale);
        }
    }

    printf("Printing evec\n");
    FILE *fh_evec = fopen(argv[4], "w");
    kjg_gsl_evec_fprintf(fh_evec, eval, evec, "%g");
}
