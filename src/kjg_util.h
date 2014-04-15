/*
 * kjg_util.h
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#ifndef KJG_UTIL_H_
#define KJG_UTIL_H_

#include <stdio.h>

/**
 * Open a file given a prefix and suffix.
 *
 * @param prefix
 * @param suffix
 * @param opentype passed to fopen
 */

FILE* kjg_fopen_suffix (const char* prefix, const char* suffix, const char* opentype);

#endif /* KJG_UTIL_H_ */
