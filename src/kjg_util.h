/*
 * kjg_util.h
 *
 *  Created on: Aug 1, 2013
 *      Author: kjg063
 */

#ifndef KJG_UTIL_H_
#define KJG_UTIL_H_

/**
 * Open a file given a prefix and suffix.
 *
 * @param prefix
 * @param suffix
 * @param opentype passed to fopen
 */

FILE* kjg_fopen_suffix (const char* prefix, const char* suffix, const char* opentype);


/**
 * Return number of seconds elapsed from two timespec structs
 *
 * @param const struct timespec t1 - start time
 * @param const struct timespec t2 - end time
 * @return double - seconds elapsed
 */

double kjg_deltat(const struct timespec t1, const struct timespec t2);

#endif /* KJG_UTIL_H_ */
