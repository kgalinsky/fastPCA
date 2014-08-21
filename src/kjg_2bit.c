#include "kjg_2bit.h"

#include <string.h>

// macros for packing array generation
#define P1(a, b, c) { KJG_2BIT_PACK(a,b,c,0), \
                      KJG_2BIT_PACK(a,b,c,1), \
                      KJG_2BIT_PACK(a,b,c,2), \
                      KJG_2BIT_PACK(a,b,c,3) }

#define P2(a, b)    { P1(a,b,0), P1(a,b,1), P1(a,b,2), P1(a,b,3) }
#define P3(a)       { P2(a,0),   P2(a,1),   P2(a,2),   P2(a,3) }
#define P4          { P3(0),     P3(1),     P3(2),     P3(3) }

// pack lookup array (faster than running PACK)
const uint8_t KJG_2BIT_PACK_LOOKUP[4][4][4][4] = P4;

// macros for unpacking array generation - instead of a for loop
#define U1(n) KJG_2BIT_UNPACK(n), \
              KJG_2BIT_UNPACK(n + 1), \
              KJG_2BIT_UNPACK(n + 2), \
              KJG_2BIT_UNPACK(n + 3)

#define U2(n) U1(n), U1(n +  4), U1(n +   8), U1(n +  12)
#define U3(n) U2(n), U2(n + 16), U2(n +  32), U2(n +  48)
#define U4(n) U3(n), U3(n + 64), U3(n + 128), U3(n + 192)
// unpack lookup array
const uint8_t KJG_2BIT_UNPACK_LOOKUP[256][4] = { U4(0) };

/**
 * Packs an array of 4 integers into one byte
 * @param *unpacked array of 4 integers
 * @return byte containing the packed integers
 */

inline uint8_t kjg_2bit_pack_unit (const uint8_t* unpacked) {
    return (KJG_2BIT_PACK_LOOKUP[unpacked[0]][unpacked[1]][unpacked[2]][unpacked[3]]);
}

size_t kjg_2bit_pack (const size_t n, const uint8_t* unpacked, uint8_t* packed) {
    size_t i, j = 0;

    // pack the whole chunks
    for (i = 0; i < n - 4; i += 4)
        packed[j++] = kjg_2bit_pack_unit(&unpacked[i]);

    // pack the last chunk
    uint8_t remainder[4] = { 0, 0, 0, 0 };
    for (; i < n; i++)
        remainder[i % 4] = unpacked[i];

    packed[j] = kjg_2bit_pack_unit(remainder);

    return (j);
}

size_t kjg_2bit_unpack (
        const size_t n,
        const uint8_t* packed,
        uint8_t* unpacked) {
    size_t i, j = 0;
    for (i = 0; i < n - 4; i += 4) {
        memcpy(&unpacked[i], &KJG_2BIT_UNPACK_LOOKUP[packed[j++]], 4);
    }
    memcpy(&unpacked[i], &KJG_2BIT_UNPACK_LOOKUP[packed[j]], n - i);
    return (j);
}

size_t kjg_2bit_unpack_or (
        const size_t n,
        const uint8_t* packed,
        const uint8_t* mask,
        uint8_t* unpacked) {
    size_t i, j = 0;
    for (i = 0; i < n; i += 4) {
        memcpy(&unpacked[i], &KJG_2BIT_UNPACK_LOOKUP[packed[j] | mask[j++]],
                sizeof(uint32_t));
    }
    return (j);
}

size_t kjg_2bit_packed_tda (const size_t n) {
    return ((n + 3) / 4);
}
