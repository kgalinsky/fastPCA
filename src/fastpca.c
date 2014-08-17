/*
 * fastpca.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_linalg.h"

#include "kjg_geno.h"
#include "kjg_genoIO.h"
#include "kjg_gsl.h"
#include "kjg_util.h"
#include "kjg_fpca.h"

//#include <lapacke.h>

// options and arguments
size_t I = 10;
size_t L = 10;
size_t K = 5;
char *GENO_FILENAME, *OUTPUT_PREFIX;

// argument parsing
void parse_args (int argc, char **argv);
extern int opterr, optopt, optind;
extern char *optarg;

int timelog (const char* message);
struct timespec t0;

struct timespec elapsed ();

int main (int argc, char **argv) {
    clock_gettime(CLOCK_REALTIME, &t0);
    char message[256];

    parse_args(argc, argv);

    FILE *fh_geno = fopen(GENO_FILENAME, "r");

    size_t n = kjg_genoIO_num_ind(fh_geno);
    size_t m = kjg_genoIO_num_snp(fh_geno, n);

    kjg_geno* X = kjg_geno_alloc(m, n);
    double* M = malloc(sizeof(double) * m);

    gsl_vector* eval = gsl_vector_alloc(K);
    gsl_matrix* evec = gsl_matrix_alloc(n, K);

    // PREP - read genotype file into memory
    sprintf(message, "Reading geno (%dx%d)", m, n);
    timelog(message);
    kjg_genoIO_fread_geno(X, fh_geno);
    fclose(fh_geno);

    // calculate the SNP means
    kjg_geno_row_means(X, M);

    timelog("fastPCA started");
    kjg_fpca(X, M, eval, evec, L, I);
    timelog("fastPCA completed");

    FILE *fh_evec = kjg_fopen_suffix(OUTPUT_PREFIX, "evec", "w");
    kjg_gsl_evec_fprintf(fh_evec, eval, evec, "%g");
    fclose(fh_evec);

    timelog("Done");
    return (0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:r:")) != -1) {
        switch (c) {
        case 'i':
            I = atoi(optarg);
            break;
        case 'k':
            K = atoi(optarg);
            break;
        case 'l':
            L = atoi(optarg);
            break;
        case 'o':
            OUTPUT_PREFIX = optarg;
            break;
        case 'r':
            KJG_FPCA_ROWS = atoi(optarg);
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

    if (K >= L) {
        fprintf(stderr, "K >= L\n");
        exit(1);
    }

    if (OUTPUT_PREFIX == NULL) {
        OUTPUT_PREFIX = GENO_FILENAME;
    }
}

int timelog (const char* message) {
    struct timespec ts = elapsed();
    return (printf("[%06d.%09d] %s\n", ts.tv_sec, ts.tv_nsec, message));
}

struct timespec elapsed () {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    if (ts.tv_nsec < t0.tv_nsec) {
        ts.tv_nsec = 1000000000 + ts.tv_nsec - t0.tv_nsec;
        ts.tv_sec--;
    }
    ts.tv_sec -= t0.tv_sec;
    return (ts);
}
