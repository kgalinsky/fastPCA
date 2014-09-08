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

    // Allocate file handles
    FILE *fh_evec = fopen(argv[1], "r");
    kjg_genoIO *gp = kjg_genoIO_fopen(argv[2], "r");
    FILE *fh_weight = fopen(argv[3], "w");

    // Get dimensions of evec file
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

    // only reading 1024 lines at a time
    size_t i, m = 1024, mr;

    // allocate vectors
    gsl_vector *eval = gsl_vector_alloc(k);
    gsl_matrix *evec = gsl_matrix_alloc(gp->n, k);
    gsl_matrix *weights = gsl_matrix_alloc(m, k);

    // ead evec
    printf("Reading evec (%dx%d)\n", gp->n, k);
    int r = kjg_gsl_evec_fscanf(fh_evec, eval, evec);
    if (r) {
        printf("Error reading evec file\n");
        exit(1);
    }
    fclose(fh_evec);

    // print eval
    fprintf(fh_weight, "#");
    fprintf(fh_weight, "%g", gsl_vector_get(eval, 0));
    for (i = 1; i < eval->size; i++) {
        fprintf(fh_weight, "\t");
        fprintf(fh_weight, "%g", gsl_vector_get(eval, i));
    }

    // calculate scale
    double *scale = malloc(sizeof(double) * eval->size);
    for (i = 0; i < k; i++) {
        scale[i] = 1.0 / sqrt(gp->m * gsl_vector_get(eval, i));
    }

    gsl_vector_free(eval);

    // do it in chunks
    kjg_geno* X = kjg_geno_alloc(1024, gp->n);
    double *M = malloc(sizeof(double) * X->m);

    size_t chunk = 0;
    while (1) {
        if (!chunk++ % 1000) fprintf(stderr, "Chunk %d\n", i);

        mr = kjg_genoIO_fread_chunk(gp, X);
        kjg_geno_row_means(X, M);

        kjg_fpca_XA(X, M, evec, weights);

        for (i = 0; i < k; i++) {
            gsl_vector_view V = gsl_matrix_column(weights, i);
            gsl_vector_scale(&V.vector, scale[i]);
        }

        if (mr == m) {
            kjg_gsl_matrix_fprintf(fh_weight, weights, "%g");
        }
        else {
            gsl_matrix_view weight_view = gsl_matrix_submatrix(weights, 0, 0,
                    mr, weights->size2);
            kjg_gsl_matrix_fprintf(fh_weight, &weight_view.matrix, "%g");
            break;
        }
    }

    kjg_geno_free(X);
    free(M);
    gsl_matrix_free(evec);
    gsl_matrix_free(weights);
}
