/*
 * kjg_util.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kjg_util.h"

FILE* kjg_fopen_suffix (const char* prefix, const char* suffix, const char* opentype) {
    char *filename = malloc(sizeof(char) * (strlen(prefix) + strlen(suffix) + 1));

    strcpy(filename, prefix);
    strcat(filename, ".");
    strcat(filename, suffix);

    FILE* fh = fopen(filename, opentype);
    free(filename);
    return(fh);
}
