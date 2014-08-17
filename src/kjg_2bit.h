/** @file kjg_2bit.h
 * @brief Packs/unpacks arrays of 2-bit integers
 * This module does the guts of packing/unpacking the 2-bit stored genotypes.
 * It is agnostic to the underlying data structure - it just wants input and
 * output arrays.
 */

#ifndef KJG_2BIT_H_
#define KJG_2BIT_H_

#include <stddef.h>
#include <stdint.h>

/** Packs an array of integers into 4-integer-per-byte array
 * @param n number of integers
 * @param *unpacked array of integers
 * @param *packed array to store packed integers
 * @return length of packed array
 */
size_t kjg_2bit_pack (const size_t n, const uint8_t* unpacked, uint8_t* packed);

/** Unpacks a 4-integer-per-byte array
 * @param n number of integers
 * @param *packed array of packed integers
 * @param *unpacked array to store unpacked integers
 * @return length of packed array
 */

size_t kjg_2bit_unpack (
        const size_t n,
        const uint8_t* packed,
        uint8_t* unpacked);
/** Unpacks a 4-integer-per-byte array with a bitwise "or" mask
 * @param n number of integers
 * @param *packed array of packed integers
 * @param *mask mask array
 * @param *unpacked array to store unpacked integers
 * @return length of packed array
 */
size_t kjg_2bit_unpack_or (
        const size_t n,
        const uint8_t* packed,
        const uint8_t* mask,
        uint8_t* unpacked);

/** Determines the length of the packed array
 * @param n length of unpacked array
 * @return Length of packed array
 */

size_t kjg_2bit_packed_tda (const size_t n);

#endif
