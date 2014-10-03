#include <stdio.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kjg_geno.h"

void write_bim (char* prefix, size_t m);
void write_fam (char* prefix, size_t n, size_t p);
void write_bed (char* prefix, size_t m, size_t n, size_t p, double FST);

void main (int argc, char** argv) {
    if (argc != 6) {
        printf("r6 <M> <N> <populations> <FST> <output prefix>\n");
        return;
    }

    size_t m = atoi(argv[1]);
    size_t n = atoi(argv[2]);
    size_t p = atoi(argv[3]);
    double FST = atof(argv[4]);
    char* prefix = argv[5];

    write_bim(prefix, m);
    write_fam(prefix, n, p);
    write_bed(prefix, m, n, p, FST);
}

void write_bim (char* prefix, size_t m) {
    char* filename = malloc(strlen(prefix) + 5);
    sprintf(filename, "%s.bim", prefix);
    FILE* fp = fopen(filename, "w");
    free(filename);

    size_t i;
    for (i = 1; i <= m; i++)
        fprintf(fp, "0\tsnp%08d\t0\t%d000\tA\tC\n", i, i);

    fclose(fp);
}

void write_fam (char* prefix, size_t n, size_t p) {
    char* filename = malloc(strlen(prefix) + 5);
    sprintf(filename, "%s.fam", prefix);
    FILE* fp = fopen(filename, "w");
    free(filename);

    size_t i, j;
    for (i = 1; i <= p; i++)
        for (j = 1; j <= n; j++)
            fprintf(fp, "pop%d\tind%05d\t0\t0\t0\t0\n", i, j);

    fclose(fp);
}

void write_bed (char* prefix, size_t m, size_t n, size_t p, double FST) {
    char* filename = malloc(strlen(prefix) + 5);
    sprintf(filename, "%s.bed", prefix);
    FILE* fp = fopen(filename, "w");
    free(filename);

    char magic[3] = { 0x6c, 0x1b, 0x01 };
    fwrite(magic, 1, 3, fp);

    gsl_rng_env_setup();

    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);

    double F = (1 - FST) / FST;

    size_t i, j, k, l;
    for (i = 0; i < m; i++) {
        double anc = gsl_rng_uniform(r)*0.9 + 0.05;
        double a = anc * F;
        double b = (1 - anc) * F;

        for (j = 0; j < p; j++) {
            double pop = gsl_ran_beta(r, a, b);

            // quicker than calling gsl_ran_binomial a million times
            double cut1 = pop * pop;
            double cut2 = 1 - (1 - pop) * (1 - pop);

            uint8_t gen[4];
            uint8_t g;

            for (k = 0; k < n; k += 4) {
                for (l = 0; l < 4; l++) {
                    double ind = gsl_rng_uniform(r);
                    gen[l] = ind < cut1 ? 0 : ind < cut2 ? 2 : 3;
                }
                g = kjg_geno_pack_unit(gen);
                fwrite(&g, 1, 1, fp);
            }
        }
    }

    fclose(fp);
}
