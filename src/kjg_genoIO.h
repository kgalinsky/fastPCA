/*
 * @file kjg_genoIO.h
 * @brief Reads geno files
 */

#ifndef KJG_GENOIO_H_
#define KJG_GENOIO_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "kjg_geno.h"

// Struct for reading geno files

typedef struct {
    const size_t m;   // number of SNPs
    const size_t n;   // number of samples
    FILE* stream;
} kjg_genoIO;

/**
 * Opens a geno file
 * @param *path path to file
 * @param *mode mode to open (supports only read for now)
 * @return point to kjg_genoIO struct
 */

kjg_genoIO* kjg_genoIO_fopen (const char* path, const char* mode);

/**
 * Closes a geno file
 * @param *gp pointer to kjg_genoIO struct
 */

int kjg_genoIO_fclose (kjg_genoIO* gp);

/**
 * Determine the number of individuals in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @return number of individuals
 */
size_t kjg_genoIO_num_ind (FILE* stream);

/**
 * Determine the number of SNPs in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @param n number of individuals
 * @return number of SNPs
 */
size_t kjg_genoIO_num_snp (FILE* stream, size_t n);

/**
 * Convert a character buffer to geno.
 *
 * @param *buffer character buffer
 * @param *x genotype array
 * @param n number of individuals
 */
void kjg_genoIO_char2int (const char* buffer, uint8_t* x, const size_t n);

/**
 * Read line from genotype file.
 *
 * @param *buffer character buffer
 * @param *x genotype array
 * @param n number of individuals
 * @param *stream genotype file object
 * @return number of characters read (should be n+1)
 */
size_t kjg_genoIO_fread (char* buffer, uint8_t* x, const size_t n, FILE* stream);

/**
 * Read geno file into struct
 *
 * @param *g pre-allocated geno struct
 * @param *stream genotype file object
 */

void kjg_genoIO_fread_geno (kjg_geno* g, FILE* stream);

#endif /* KJG_GENOIO_H_ */
