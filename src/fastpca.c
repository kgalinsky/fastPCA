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
#include "kjg_bedIO.h"
#include "kjg_gsl.h"
#include "kjg_util.h"
#include "kjg_fpca.h"

//#include <lapacke.h>

// options and arguments
size_t I = 10;
size_t L = 20;
size_t K = 10;
int B = 0;
char *GENO_FILENAME, *OUTPUT_PREFIX;
char *SNP_MASK, *IND_MASK;

// argument parsing
void parse_args (int argc, char **argv);
void print_usage (const char *message);

extern int opterr, optopt, optind;
extern char *optarg;

int timelog (const char* message);
struct timespec t0;
struct timespec elapsed ();

uint8_t *read_mask (size_t n, const char *filename);

int main (int argc, char **argv) {
    clock_gettime(CLOCK_REALTIME, &t0);
    char message[256];

    parse_args(argc, argv);

    // PREP - read genotype file into memory
    kjg_geno* X;
    if (B) {
        kjg_bedIO *bp = kjg_bedIO_bfile_fopen(GENO_FILENAME, "r");
        if (!bp) print_usage("Unable to open bed/bim/fam");

        uint8_t* SNPmask = 0;
        if (SNP_MASK) {
            SNPmask = read_mask(bp->m, SNP_MASK);
            if (!SNPmask) print_usage("Unable to read SNP mask");
        }

        uint8_t* indmask = 0;
        if (IND_MASK) {
            indmask = read_mask(bp->m, IND_MASK);
            if (!indmask) print_usage("Unable to read individual mask");
        }

        sprintf(message, "Reading geno (%dx%d)", bp->m, bp->n);
        timelog(message);

        X = kjg_bedIO_fread_geno(bp, SNPmask, indmask);
        kjg_bedIO_fclose(bp);

        if (SNPmask) free(SNPmask);
        if (indmask) free(indmask);

        sprintf(message, "Finished reading geno (%dx%d)", X->m, X->n);
        timelog(message);
    }
    else {
        kjg_genoIO *gp = kjg_genoIO_fopen(GENO_FILENAME, "r");
        if (!gp) {
            print_usage("Unable to open geno");
        }

        sprintf(message, "Reading geno (%dx%d)", gp->m, gp->n);
        timelog(message);

        X = kjg_genoIO_fread_geno(gp);
        kjg_genoIO_fclose(gp);
    }

    // calculate the SNP means
    timelog("Calculating SNP allele frequencies");
    kjg_geno_set_norm(X, 0);

    // run fast PCA
    gsl_vector* eval = gsl_vector_alloc(K);
    gsl_matrix* evec = gsl_matrix_alloc(X->n, K);

    timelog("fastPCA started");
    kjg_fpca(X, eval, evec, L, I);
    timelog("fastPCA completed");

    FILE *fh_evec = kjg_util_fopen_suffix(OUTPUT_PREFIX, "evec", "w");
    kjg_gsl_evec_fprintf(fh_evec, eval, evec, "%g");
    fclose(fh_evec);

    timelog("Done");
    return (0);
}

void parse_args (int argc, char **argv) {
    int c;

    opterr = 1;
    while ((c = getopt(argc, argv, "i:k:l:o:bm:n:")) != -1) {
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
        case 'b':
            B = 1;
            break;
        case 'm':
            SNP_MASK = optarg;
            break;
        case 'n':
            IND_MASK = optarg;
            break;
        default:
            print_usage(0);
        }
    }

    if (optind == argc - 1) {
        GENO_FILENAME = argv[optind];
    }
    else if (optind == argc) print_usage("fastpca: No input file specified");
    else print_usage("fastpca: Too many arguments");

    if (K >= L) print_usage("fastpca: K >= L");

    if (OUTPUT_PREFIX == NULL) {
        OUTPUT_PREFIX = GENO_FILENAME;
    }
}

void print_usage (const char *message) {
    if (message) fprintf(stderr, "fastpca: %s\n", message);
    fprintf(stderr,
            "Usage: fastpca [options] <input>\n"
                    "Options:\n"
                    "  -k <PCs>             Number of PCs to compute (K). Default is 10.\n"
                    "  -l <width>           Width of random matrix (L). L > K. Default is 20.\n"
                    "  -i <iterations>      Number of iterations (I). Default is 10.\n"
                    "  -o <output prefix>   Output prefix. Default is to use input.\n"
                    "  -b                   Tells fastpca that we are using a bed file.\n"
                    "  -m <SNP mask>        File with SNP masking (1 means don't use).\n"
                    "  -n <ind mask>        File with individual masking (1 means don't use).\n"
                    "  <input>              Input geno file or bed/bim/fam prefix with -b.\n");
    exit(1);
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

uint8_t *read_mask (size_t n, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return (NULL);

    uint8_t *mask = calloc(n, 1);
    size_t i;
    for (i = 0; i < n; i++)
        if (fscanf(fp, "%d", &mask[i]) == EOF) {
            free(mask);
            fclose(fp);
            return (NULL);
        }
    fclose(fp);

    return (mask);
}

