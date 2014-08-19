#include "kjg_genoIO.h"

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

// Map characters to integer
static const uint8_t KJG_GENOIO_CHAR_MAP[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        // The magic happens here
        0, 1, 2, 4, 4, 4, 4, 4, 4, 3,
        // Maps characters to the appropriate numbers
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4 };

kjg_genoIO* kjg_genoIO_fopen (const char* path, const char* mode) {
    if (mode[0] != 'r') return (NULL); // TODO support writing

    FILE* stream = fopen(path, mode);
    if (stream == NULL) return (NULL);

    size_t n = kjg_genoIO_num_ind(stream);
    size_t m = kjg_genoIO_num_snp(stream, n);

    kjg_genoIO pre = { m, n, stream };
    kjg_genoIO* gp = malloc(sizeof(kjg_genoIO));
    memcpy(gp, &pre, sizeof(kjg_genoIO));

    return (gp);
}

int kjg_genoIO_fclose (kjg_genoIO* gp) {
    int r = fclose(gp->stream);
    free(gp);
    return (r);
}

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
        kjg_geno_set_row(g, i, x);
    }

    free(buffer);
    free(x);
}
