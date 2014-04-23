/*
 * kjg_util.c
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

double kjg_deltat(const struct timespec t1, const struct timespec t2) {
    return(t2.tv_sec - t1.tv_sec + (double)(t2.tv_nsec - t1.tv_nsec)/(1000000000));
}
