/**
 * @file kjg_2bit.h
 * @brief Packs/unpacks arrays of 2-bit integers
 * This module does the guts of packing/unpacking the 2-bit stored genotypes.
 * It is agnostic to the underlying data structure - it just wants input and
 * output arrays.
 */

#include <stddef.h>
#include <stdint.h>

/**
 * Unpacks a 4-integer-per-byte array
 * @param n number of integers
 * @param *packed array of packed integers
 * @param *unpacked array to store unpacked integers
 * @return length of packed array
 */

size_t kjg_2bit_unpack(const size_t n, const uint8_t* packed, uint8_t* unpacked);
/**
 * Unpacks a 4-integer-per-byte array with a bitwise "or" mask
 * @param n number of integers
 * @param *packed array of packed integers
 * @param *mask mask array
 * @param *unpacked array to store unpacked integers
 * @return length of packed array
 */
size_t kjg_2bit_unpack_or(const size_t n, const uint8_t* packed,
		const uint8_t* mask, uint8_t* unpacked);

/**
 * Packs an array of integers into 4-integer-per-byte array
 * @param n number of integers
 * @param *unpacked array of integers
 * @param *packed array to store packed integers
 * @return length of packed array
 */
size_t kjg_2bit_pack(const size_t n, const uint8_t* unpacked, uint8_t* packed);

/**
 * Packs an array of 4 integers into one byte
 * @param *unpacked array of 4 integers
 * @return Byte containing the packed integers
 */
uint8_t kjg_2bit_pack_unit(const uint8_t* unpacked) {
	return (PACK_LOOKUP[unpacked[3]][unpacked[2]][unpacked[1]][unpacked[0]]);
}

/**
 * Determines the length of the packed array
 * @param n length of unpacked array
 * @return Length of packed array
 */

size_t kjg_2bit_packed_tda(const size_t n) {
	return ((n + 3) / 4);
}
