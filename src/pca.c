/*
 * pca.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <time.h>

#include <lapacke.h>

#include "kjg_geno.h"
#include "kjg_genoIO.h"
#include "kjg_GRM.h"
#include "kjg_util.h"

// options and arguments
size_t K = 10;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args (int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

// helper functions
void scale_evals (double* evals, const size_t n);
void fprintf_evals (
        FILE* fh_eval,
        const char* format,
        const double* evals,
        const size_t n);
void fprintf_evecs (
        FILE* fh_evec,
        const char* format,
        const double* evals,
        const double* evecs,
        const size_t n,
        const size_t K);

// MAIN FUNCTION
int main (int argc, char** argv) {
    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");
//    FILE *fh_eval = kjg_fopen_suffix(OUTPUT_PREFIX, "eval", "w");
    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    double* GRM = kjg_GRM_init(n);

    {
        char* buffer = malloc(sizeof(char) * (n + 1));
        uint8_t* x = malloc(sizeof(uint8_t) * n);

        size_t i;
        for (i = 0; i < m; i++) {
            kjg_genoIO_fread(buffer, x, n, fh_geno);
            kjg_GRM_update(x, GRM, n);
        }

        free(buffer);
        free(x);
    }

    double* evals = malloc(sizeof(double) * n);   // eigenvalues
    double* evecs = malloc(sizeof(double) * n * n); // eigenvectors

    LAPACKE_dspevd(LAPACK_COL_MAJOR, 'V', 'L', n, GRM, evals, evecs, n);

    free(GRM);

    scale_evals(evals, n);
//    fprintf_evals(fh_eval, "%g", evals, n);
    fprintf_evecs(fh_evec, "%g", evals, evecs, n, K);

    free(evals);
    free(evecs);

    return (0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:")) != -1) {
        switch (c) {
        case 'k':
            K = atoi(optarg);
            break;
        case 'o':
            OUTPUT_PREFIX = optarg;
            break;
        default:
            fprintf(stderr, "Unrecognized option '-%c'\n", optopt);
            exit(1);
        }
    }

    if (optind == argc - 1) {
        GENO_FILENAME = argv[optind];
    }
    else if (optind == argc) {
        fprintf(stderr, "No input file specified\n");
        exit(1);
    }
    else {
        fprintf(stderr, "Too many arguments\n");
        exit(1);
    }

    if (OUTPUT_PREFIX == NULL) {
        OUTPUT_PREFIX = GENO_FILENAME;
    }
}

void scale_evals (double* evals, const size_t n) {
    size_t i;
    double sum = 0;
    for (i = 0; i < n; i++) {
        sum += evals[i];
    }
    double scale = (double) n / sum;
    for (i = 0; i < n; i++) {
        evals[i] *= scale;
    }
}

void fprintf_evals (
        FILE* fh_eval,
        const char* format,
        const double* evals,
        const size_t n) {
    size_t i;
    for (i = n; i-- > 0;) {
        fprintf(fh_eval, format, evals[i]);
        fprintf(fh_eval, "\n");
    }
}

void fprintf_evecs (
        FILE* fh_evec,
        const char* format,
        const double* evals,
        const double* evecs,
        const size_t n,
        const size_t K) {
    size_t i, j;
    fprintf(fh_evec, "#");
    fprintf(fh_evec, format, evals[n - 1]);
    for (i = n - 1; i-- > n - K;) {
        fprintf(fh_evec, "\t");
        fprintf(fh_evec, format, evals[i]);
    }
    fprintf(fh_evec, "\n");
    for (i = 0; i < n; i++) {
        fprintf(fh_evec, format, evecs[n * (n - 1) + i]);
        for (j = n - 2; j-- > n - K;) {
            fprintf(fh_evec, "\t");
            fprintf(fh_evec, format, evecs[n * j + i]);
        }
        fprintf(fh_evec, "\n");
    }
}
