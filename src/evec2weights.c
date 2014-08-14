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
    FILE *fh_geno = fopen(argv[2], "r");

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

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    gsl_vector *eval = gsl_vector_alloc(k);
    gsl_matrix *evec = gsl_matrix_alloc(n, k);
    gsl_matrix *weights = gsl_matrix_alloc(m, k);

    kjg_geno *X = kjg_geno_alloc(m, n);
    double *M = malloc(sizeof(double) * m);

    printf("Reading evec (%dx%d)\n", n, k);
    int r = kjg_gsl_evec_fscanf(fh_evec, eval, evec);
    if (r) {
        printf("Error reading evec file\n");
        exit(1);
    }
    fclose(fh_evec);

    printf("Reading geno (%dx%d)\n", m, n);
    kjg_genoIO_fread_geno(X, fh_geno);
    fclose(fh_geno);

    printf("Computing weights (%dx%d)\n", m, k);
    kjg_geno_row_means(X, M);
    kjg_fpca_XG(X, M, evec, weights);

    printf("Scaling weights\n");
    {
        size_t i;
        for (i = 0; i < k; i++) {
            double scale = 1.0 / sqrt(m * gsl_vector_get(eval, i));
            gsl_vector_view V = gsl_matrix_column(weights, i);
            gsl_vector_scale(&V.vector, scale);
        }
    }

    printf("Printing weights\n");
    FILE *fh_weight = fopen(argv[3], "w");
    kjg_gsl_evec_fprintf(fh_weight, eval, weights, "%g");
}
