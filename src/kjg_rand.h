/*
 * kjg_rand.h
 *
 *  Created on: Apr 23, 2014
 *      Author: Kevin
 */

#ifndef KJG_RAND_H_
#define KJG_RAND_H_

/**
 * Generate a single uniform(0,1) random number.
 */

double kjg_runif();

/**
 * Generate n normal(0,1) random numbers.
 * @param const size_t n - number of normals to generate
 * @param double* Z      - array to store normal
 */

void kjg_rnorms(const size_t n, double* Z);

#endif /* KJG_RAND_H_ */
