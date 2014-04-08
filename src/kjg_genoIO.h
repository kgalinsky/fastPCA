/*
 * kjg_genoIO.h
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 */

#ifndef KJG_GENOIO_H_
#define KJG_GENOIO_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "kjg_geno.h"

/**
 * Determine the number of individuals in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @return number of individuals
 */
size_t kjg_genoIO_num_ind(FILE* stream);

/**
 * Determine the number of SNPs in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @param n number of individuals
 * @return number of SNPs
 */
size_t kjg_genoIO_num_snp(FILE* stream, size_t n);

/**
 * Convert a character buffer to geno.
 *
 * @param *buffer character buffer
 * @param *x genotype array
 * @param n number of individuals
 */
void kjg_genoIO_char2int(const char* buffer, uint8_t* x, const size_t n);

/**
 * Read line from genotype file.
 *
 * @param *buffer character buffer
 * @param *x genotype array
 * @param n number of individuals
 * @param *stream genotype file object
 * @return number of characters read (should be n+1)
 */
size_t kjg_genoIO_fread(char* buffer, uint8_t* x, const size_t n, FILE* stream);

/**
 * Read geno file into struct
 *
 * @param *buffer character buffer to store lines
 * @param *x integer buffer to store unpacked genotypes
 * @param *g pre-allocated geno struct
 * @param *stream genotype file object
 */

void kjg_genoIO_fread_geno(char* buffer, uint8_t* x, kjg_geno* g, FILE* stream);

#endif /* KJG_GENOIO_H_ */
