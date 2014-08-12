/**
 * @file kjg_2bit.h
 * @brief Packs/unpacks arrays of 2-bit integers
 */

#include <stddef.h>
#include <stdint.h>

/**
 * Packs an array of integers into 4-integer-per-byte array
 * @param n number of integers
 * @param *unpacked array of integers
 * @param *packed array to store packed integers
 */

void kjg_2bit_pack(const size_t n, const uint8_t* unpacked, uint8_t* packed);

/**
 * Unpacks a 4-integer-per-byte array
 * @param n number of integers
 * @param *packed array of packed
 * @param *unpacked array to store unpacked integers
 */

void kjg_2bit_unpack(const size_t n, const uint8_t* packed, uint8_t* unpacked);
