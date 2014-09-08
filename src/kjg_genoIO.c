#include "kjg_genoIO.h"

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

size_t KJG_GENOIO_MAX_BUFFER = 1048576; // 1MB buffer for reading

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

kjg_geno* kjg_genoIO_fread_geno (kjg_genoIO* gp) {
    kjg_geno* g = kjg_geno_alloc(gp->m, gp->n);

    size_t n1 = gp->n + 1;                  // n + 1
    size_t nb = KJG_GENOIO_MAX_BUFFER / n1; // number of lines in the buffer

    char *buffer = malloc(sizeof(char) * n1 * nb);
    uint8_t *x = malloc(sizeof(uint8_t) * gp->n);

    size_t i, j = 0;
    size_t nr;                           // number of lines read
    do {
        nr = fread(buffer, sizeof(char) * n1, nb, gp->stream);
        for (i = 0; i < nr; i++) {
            kjg_genoIO_char2int(&buffer[i * n1], x, gp->n);
            kjg_geno_set_row(g, j++, x);
        }
    } while (nr == nb);

    free(buffer);
    free(x);

    return (g);
}

size_t kjg_genoIO_fread_chunk (kjg_genoIO* gp, kjg_geno* g) {
    size_t i, n1 = g->n + 1;

    char *buffer = malloc(sizeof(char) * n1 * g->m);
    uint8_t *x = malloc(sizeof(uint8_t) * g->n);

    size_t nr = fread(buffer, sizeof(char) * n1, g->m, gp->stream);

    for (i = 0; i < nr; i++) {
        kjg_genoIO_char2int(&buffer[i * n1], x, gp->n);
        kjg_geno_set_row(g, i, x);
    }

    return (nr);
}
