#include "kjg_bedIO.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

kjg_bedIO* kjg_bedIO_fopen (
        const char* path,
        const char* mode,
        const size_t m,
        const size_t n) {
    if (mode[0] != 'r') return (NULL); // TODO support writing

    FILE* stream = fopen(path, mode);
    if (stream == NULL) return (NULL);

    kjg_bedIO pre = { m, n, stream };
    kjg_bedIO* bp = malloc(sizeof(kjg_bedIO));
    memcpy(bp, &pre, sizeof(kjg_bedIO));

    char magic[3];
    fread(magic, 1, 3, stream);
    if ((magic[0] != 0x6c) || (magic[1] != 0x1b) || (magic[2] != 0x01)) {
        fprintf(stderr, "Bad magic numbers: %02x %02x %02x\n", magic[0],
                magic[1], magic[2]);
        exit(1);
    }
    return (bp);
}

int kjg_bedIO_fclose (kjg_bedIO* bp) {
    int r = fclose(bp->stream);
    free(bp);
    return (r);
}

kjg_geno* kjg_bedIO_fread_geno (kjg_bedIO* bp) {
    kjg_geno* g = kjg_geno_alloc(bp->m, bp->n);
    fread(g->data, 1, g->tda * g->m, bp->stream);
    return (g);
}
