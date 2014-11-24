#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kjg_geno.h"

void
parse_args (int argc, char** argv);

void
write_bim (char* prefix, size_t m);
void
write_fam (char* prefix, size_t n, size_t p);
void
write_bed (char* prefix, size_t m, size_t n, size_t p, double FST, double *MAF);

size_t M = 10000, N = 1000, P = 2;
double FST = 0.01, MAF[2] =
  { 0.05, 0.5 };
char* PFX = NULL;

void
main (int argc, char** argv)
{
  parse_args (argc, argv);

  write_bim (PFX, M);
  write_fam (PFX, N, P);
  write_bed (PFX, M, N, P, FST, MAF);
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
      { 0, 0, 0, 0 } };

  int help_flag = 0, c, option_index = 0;

  while ((c = getopt_long (argc, argv, "m:n:f:a:p:", long_options,
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
          M = strtoul (optarg, 0, 0);
          break;

          // Samples
        case 'n':
          sscanf (optarg, "%ux%u", &N, &P);
          if (P == 0)
            P = 1;
          break;

          //FST
        case 'f':
          FST = strtof(optarg, 0, 0);
          break;

          //MAF
        case 'a':
          sscanf (optarg, "%g-%g", &MAF[0], &MAF[1]);
          if ((MAF[0] < 0) || (MAF[0] > MAF[1]) || (MAF[1] > 0.5))
            {
              fprintf (stderr, "MAF bad: %g-%g\n", MAF[0], MAF[1]);
              help_flag = 1;
            }
          break;

          //PREFIX
        case 'p':
          asprintf (&PFX, "%s", optarg);
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

  if (PFX == NULL)
    asprintf (&PFX, "star.%dx%dx%d_%g_%g-%g", M, N, P, FST, MAF[0], MAF[1]);
}

void
write_bim (char* prefix, size_t m)
{
  char* filename = malloc (strlen (prefix) + 5);
  sprintf (filename, "%s.bim", prefix);
  FILE* fp = fopen (filename, "w");
  free (filename);

  size_t i;
  for (i = 1; i <= m; i++)
    fprintf (fp, "0\tsnp%08d\t0\t%d000\tA\tC\n", i, i);

  fclose (fp);
}

void
write_fam (char* prefix, size_t n, size_t p)
{
  char* filename = malloc (strlen (prefix) + 5);
  sprintf (filename, "%s.fam", prefix);
  FILE* fp = fopen (filename, "w");
  free (filename);

  size_t i, j;
  for (i = 1; i <= p; i++)
    for (j = 1; j <= n; j++)
      fprintf (fp, "pop%d\tind%05d\t0\t0\t0\t0\n", i, j);

  fclose (fp);
}

void
write_bed (char* prefix, size_t m, size_t n, size_t p, double FST, double *MAF)
{
  char* filename = malloc (strlen (prefix) + 5);
  sprintf (filename, "%s.bed", prefix);
  FILE* fp = fopen (filename, "w");
  free (filename);

  char magic[3] =
    { 0x6c, 0x1b, 0x01 };
  fwrite (magic, 1, 3, fp);

  gsl_rng_env_setup ();

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc (T);

  double F = (1 - FST) / FST;

  size_t i, j, k, l;
  for (i = 0; i < m; i++)
    {
      double anc = gsl_ran_flat (r, MAF[0], MAF[1]);
      if (gsl_ran_bernoulli (r, 0.5))
        anc = 1 - anc;

      double a = anc * F;
      double b = (1 - anc) * F;

      for (j = 0; j < p; j++)
        {
          double pop = gsl_ran_beta (r, a, b);

          // quicker than calling gsl_ran_binomial a million times
          double cut1 = pop * pop;
          double cut2 = 1 - (1 - pop) * (1 - pop);

          uint8_t gen[4];
          uint8_t g;

          for (k = 0; k < n; k += 4)
            {
              for (l = 0; l < 4; l++)
                {
                  double ind = gsl_rng_uniform (r);
                  gen[l] = ind < cut1 ? 0 : ind < cut2 ? 2 : 3;
                }
              g = kjg_geno_pack_unit (gen);
              fwrite (&g, 1, 1, fp);
            }
        }
    }

  fclose (fp);
}
