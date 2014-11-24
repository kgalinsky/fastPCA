#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kjg_geno.h"
#include "kjg_geno_rand.h"
#include "kjg_util.h"

void
parse_args (int argc, char** argv);

void
write_bim (char* prefix, size_t m);
void
write_fam (char* prefix, size_t n, size_t p);
void
write_bed (char* prefix, size_t M, size_t N, size_t P, double FST, double *MAF);
void
write_bed_pair (char* prefix, size_t M, size_t N, size_t P, double FST,
                double *MAF, double *cor);

size_t opt_M = 10000, opt_N = 1000, opt_P = 2;
double opt_FST = 0.01, opt_MAF[2] =
  { 0.05, 0.5 }, *opt_cor = NULL;
char* opt_PFX = NULL;

void
main (int argc, char** argv)
{
  parse_args (argc, argv);

  write_bim (opt_PFX, opt_M);
  write_fam (opt_PFX, opt_N, opt_P);
  if (opt_cor)
    write_bed_pair (opt_PFX, opt_M, opt_N, opt_P, opt_FST, opt_MAF, opt_cor);
  else
    write_bed (opt_PFX, opt_M, opt_N, opt_P, opt_FST, opt_MAF);
}

void
parse_args (int argc, char** argv)
{
  static struct option long_options[] =
    {
      { "help", no_argument, 0, 'h' },
      { "M", required_argument, 0, 'm' },
      { "N", required_argument, 0, 'n' },
      { "FST", required_argument, 0, 'f' },
      { "MAF", required_argument, 0, 'a' },
      { "prefix", required_argument, 0, 'p' },
      { "cor", required_argument, 0, 'r' },
      { 0, 0, 0, 0 } };

  int help_flag = 0, c, option_index = 0;

  while ((c = getopt_long (argc, argv, "m:n:f:a:p:r:", long_options,
                           &option_index)) != -1)
    {
      switch (c)
        {
        case 0:
          break;
        case 'h':
          help_flag = 1;
          break;

          // SNPs
        case 'm':
          errno = 0;
          opt_M = strtoul (optarg, 0, 0);
          break;

          // Samples
        case 'n':
          sscanf (optarg, "%ux%u", &opt_N, &opt_P);
          if (opt_P == 0)
            opt_P = 1;
          break;

          //FST
        case 'f':
          opt_FST = strtof (optarg, 0);
          break;

          //MAF
        case 'a':
          sscanf (optarg, "%lg-%lg", &opt_MAF[0], &opt_MAF[1]);
          if ((opt_MAF[0] < 0) || (opt_MAF[0] > opt_MAF[1])
              || (opt_MAF[1] > 0.5))
            {
              fprintf (stderr, "MAF bad: %g-%g\n", opt_MAF[0], opt_MAF[1]);
              help_flag = 1;
            }
          break;

          //PREFIX
        case 'p':
          asprintf (&opt_PFX, "%s", optarg);
          break;

          //Correlation
        case 'r':
          opt_cor = malloc (2 * sizeof(double));
          int err = sscanf (optarg, "%lg,%lg", &opt_cor[0], &opt_cor[1]);
          break;
        case '?':
          break;
        default:
          abort ();
        }
    }

  if (help_flag)
    {
      fprintf (
          stderr,
          "sim_star - simulate star-shaped ancestry\n"
          "\n"
          "Usage:\n"
          "  sim_start [options]\n"
          "  sim_star -m <SNPs> -n <samples>x<populations> -f <FST>\n"
          "           -a <lower>-<upper> -p <prefix>\n"
          "\n"
          "Options:\n"
          "  -m/--M=<SNPs>            Number of SNPs (default=10k)\n"
          "  -n/--N=<samples>x<pops>  Number of samples/populations (default=1000x2)\n"
          "  -f/--FST=<FST>           FST (default=0.01)\n"
          "  -a/--MAF=<lower>-<upper> Lower/upper bounds on MAF (between 0 and 0.5\n"
          "                           - default=0.05-0.5)\n"
          "  -p/--prefix=<prefix>     Output prefix (default=star.MxNxP_FST_lower-upper\n");
      exit (1);
    }

  if (opt_PFX == NULL)
    if (opt_cor)
      asprintf (&opt_PFX, "star.%dx%dx%d_%g_%g-%g_%g,%g", opt_M, opt_N, opt_P,
                opt_FST, opt_MAF[0], opt_MAF[1], opt_cor[0], opt_cor[1]);
    else
      asprintf (&opt_PFX, "star.%dx%dx%d_%g_%g-%g", opt_M, opt_N, opt_P,
                opt_FST, opt_MAF[0], opt_MAF[1]);
}

void
write_bim (char* prefix, size_t m)
{
  FILE* fp = kjg_util_fopen_suffix (prefix, "bim", "w");

  size_t i;
  for (i = 1; i <= m; i++)
    fprintf (fp, "0\tsnp%08d\t0\t%d000\tA\tC\n", i, i);

  fclose (fp);
}

void
write_fam (char* prefix, size_t n, size_t p)
{
  FILE* fp = kjg_util_fopen_suffix (prefix, "fam", "w");

  size_t i, j;
  for (i = 1; i <= p; i++)
    for (j = 1; j <= n; j++)
      fprintf (fp, "pop%d\tind%05d\t0\t0\t0\t0\n", i, j);

  fclose (fp);
}

void
write_bed (char* prefix, size_t M, size_t N, size_t P, double FST, double *MAF)
{
  FILE* fp = kjg_util_fopen_suffix (prefix, "bed", "w");
  FILE* ld = kjg_util_fopen_suffix (prefix, "ld", "w");
  char magic[3] =
    { 0x6c, 0x1b, 0x01 };
  fwrite (magic, 1, 3, fp);

  gsl_rng_env_setup ();

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc (T);

  double F = (1 - FST) / FST;

  size_t n = P * N;
  size_t tda = (n + 3) / 4;
  double* AF = malloc (P * sizeof(double));
  uint8_t* x = malloc (n * sizeof(uint8_t));
  uint8_t* p = malloc (tda * sizeof(uint8_t));

  size_t m, j, k, l;
  for (m = 0; m < M; m++)
    {
      double anc = kjg_geno_rand_anc (r, MAF);
      kjg_geno_rand_star_AF (r, AF, anc, F, P);
      kjg_geno_rand_row (r, x, N, P, AF);
      kjg_geno_pack (n, x, p);
      fwrite (p, 1, tda, fp);
    }

  free (AF);
  free (x);

  fclose (fp);
}

void
write_bed_pair (char* prefix, size_t M, size_t N, size_t P, double FST,
                double *MAF, double *cor)
{
  FILE* fp = kjg_util_fopen_suffix (prefix, "bed", "w");
  FILE* ld = kjg_util_fopen_suffix (prefix, "ldbed", "w");

  char magic[3] =
    { 0x6c, 0x1b, 0x01 };
  fwrite (magic, 1, 3, fp);
  fwrite (magic, 1, 3, ld);

  gsl_rng_env_setup ();

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc (T);

  double F = (1 - FST) / FST;

  size_t n = P * N;
  size_t tda = (n + 3) / 4;
  double* AF = malloc (P * sizeof(double));
  uint8_t* x = malloc (n * sizeof(uint8_t));
  uint8_t* y = malloc (n * sizeof(uint8_t));
  uint8_t* p = malloc (tda * sizeof(uint8_t));

  size_t m, j, k, l;
  for (m = 0; m < M; m++)
    {
      double anc = kjg_geno_rand_anc (r, MAF);
      kjg_geno_rand_star_AF (r, AF, anc, F, P);
      kjg_geno_rand_row (r, x, N, P, AF);
      kjg_geno_rand_ld_row (r, y, N, P, AF, x, cor);
      kjg_geno_pack (n, x, p);
      fwrite (p, 1, tda, fp);
      kjg_geno_pack (n, y, p);
      fwrite (p, 1, tda, ld);
    }

  free (AF);
  free (x);
  free (y);

  fclose (fp);
  fclose (ld);
}
