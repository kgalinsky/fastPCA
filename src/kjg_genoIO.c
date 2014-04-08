/*
 * kjg_genoIO.c
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 */

#include "kjg_genoIO.h"
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

const uint8_t KJG_GENOIO_CHAR_MAP[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    0, 1, 2, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

size_t kjg_genoIO_num_ind(FILE* stream) {
	fseek(stream, 0, 0);

	size_t n = 0;
	char c;

	while (1) {
		fscanf(stream, "%c", &c);
		if (c == '\n')
			break;
		n++;
	}
	return (n);
}

size_t kjg_genoIO_num_snp(FILE* stream, size_t n) {
	fseek(stream, 0, SEEK_END);
	long int l = ftell(stream);
	fseek(stream, 0, 0);
	return (l / (n + 1));
}

void kjg_genoIO_char2int(const char* buffer, uint8_t* x, const size_t n) {
	size_t i;
	for (i = 0; i < n; i++) {
		x[i] = KJG_GENOIO_CHAR_MAP[(size_t) buffer[i]];
	}
}

size_t kjg_genoIO_fread(char* buffer, uint8_t* x, const size_t n, FILE* stream) {
	size_t r = fread(buffer, 1, n + 1, stream);
	kjg_genoIO_char2int(buffer, x, r - 1);
	return (r);
}

void kjg_genoIO_fread_geno(char* buffer, uint8_t* x, kjg_geno* g, FILE* stream) {
	size_t i;
	for (i = 0; i < g->m; i++) {
		kjg_genoIO_fread(buffer, x, g->n, stream);
		kjg_geno_set_row(x, g, i);
	}
}
