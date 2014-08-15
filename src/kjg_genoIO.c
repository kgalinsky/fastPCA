#include "kjg_genoIO.h"

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

// Map characters to integer
static const uint8_t KJG_GENOIO_CHAR_MAP[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 1, 2, 4, 4, 4, 4, 4, 4, 3,
        4, 4, 4, 4, 4,
        4, // the magic happens here
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

size_t kjg_genoIO_num_ind (FILE* stream) {
    fseek(stream, 0, 0);

    size_t n = 0;
    char c;

    while (1) {
        fscanf(stream, "%c", &c);
        if (c == '\n') break;
        n++;
    }
    return (n);
}

size_t kjg_genoIO_num_snp (FILE* stream, size_t n) {
    fseek(stream, 0, SEEK_END);
    long int l = ftell(stream);
    fseek(stream, 0, 0);
    return (l / (n + 1));
}

void kjg_genoIO_char2int (const char* buffer, uint8_t* x, const size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        x[i] = KJG_GENOIO_CHAR_MAP[(size_t) buffer[i]];
    }
}

size_t kjg_genoIO_fread (char* buffer, uint8_t* x, const size_t n, FILE* stream) {
    size_t r = fread(buffer, 1, n + 1, stream);
    kjg_genoIO_char2int(buffer, x, r - 1);
    return (r);
}

void kjg_genoIO_fread_geno (kjg_geno* g, FILE* stream) {
    char *buffer = malloc(sizeof(char) * (g->n + 1));
    uint8_t *x = malloc(sizeof(uint8_t) * g->n);

    size_t i;
    for (i = 0; i < g->m; i++) {
        kjg_genoIO_fread(buffer, x, g->n, stream);
        kjg_geno_set_row(x, g, i);
    }

    free(buffer);
    free(x);
}
