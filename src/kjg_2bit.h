/**
 * @file kjg_2bit.h
 * @brief Packs/unpacks arrays of 2-bit integers
 */

#include <stddef.h>
#include <stdint.h>

#define U1(n)     n,            n+1,              n+2,              n+3
#define U2(n) U1(n),      U1(n+256),      U1(n+256*2),      U1(n+256*3)
#define U3(n) U2(n),    U2(n+65536),    U2(n+65536*2),    U2(n+65536*3)

/**
 * Array that contains unpack lookups.
 */

static const uint32_t UNPACK_LOOKUP[256] = {
        U3(0), U3(16777216), U3(16777216*2), U3(16777216*3)
};

/**
 * Packs 4 integers into a byte.
 */
#define PACK(a, b, c, d) ((a*4 + b)*4 + c)*4 + d

#define P1(a, b, c)    { PACK(a,b,c,0), PACK(a,b,c,1), PACK(a,b,c,2), PACK(a,b,c,3) }
#define P2(a, b)       {     P1(a,b,0),     P1(a,b,1),     P1(a,b,2),     P1(a,b,3) }
#define P3(a)          {       P2(a,0),       P2(a,1),       P2(a,2),       P2(a,3) }

/**
 * Array that contains pack lookups.
 */

static const uint8_t PACK_LOOKUP[4][4][4][4] = { P3(0), P3(1), P3(2), P3(3) };

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
