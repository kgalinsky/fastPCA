#include "kjg_2bit.h"

#include <string.h>

// unpack
#define UNPACK(n) { (n) & 3, ((n) >> 2) & 3, ((n) >> 4) & 3, ((n) >> 6) & 3 }

// macros for unpacking array generation
#define U1(n) UNPACK(n),    UNPACK(n | 1),    UNPACK(n | 2),    UNPACK(n | 3)
#define U2(n)     U1(n), U1(n | (1 << 2)), U1(n | (2 << 2)), U1(n | (3 << 2))
#define U3(n)     U2(n), U2(n | (1 << 4)), U2(n | (2 << 4)), U2(n | (3 << 4))

// unpack lookup array
const uint8_t KJG_2BIT_UNPACK_LOOKUP[256][4] =
    { U3(0), U3((1 << 6)), U3((2 << 6)), U3((3 << 6)) };

// packs 4 integers into a byte
#define PACK(a, b, c, d) (d << 6) | (c << 4) | (b << 2) | a

// macros for pack array generation
#define P1(b, c, d) { PACK(0,b,c,d), PACK(1,b,c,d), PACK(2,b,c,d), PACK(3,b,c,d) }
#define P2(c, d)    {     P1(0,c,d),     P1(1,c,d),     P1(2,c,d),     P1(3,c,d) }
#define P3(d)       {       P2(0,d),       P2(1,d),       P2(2,d),       P2(3,d) }

// pack lookup array (faster than running PACK)
const uint8_t KJG_2BIT_PACK_LOOKUP[4][4][4][4] =
    { P3(0), P3(1), P3(2), P3(3) };

/** Packs an array of 4 integers into one byte
 * @param *unpacked array of 4 integers
 * @return byte containing the packed integers
 */

inline uint8_t kjg_2bit_pack_unit (const uint8_t* unpacked) {
    return (KJG_2BIT_PACK_LOOKUP[unpacked[3]][unpacked[2]][unpacked[1]][unpacked[0]]);
}

size_t kjg_2bit_pack (const size_t n, const uint8_t* unpacked, uint8_t* packed) {
    size_t i, j = 0;

    // pack the whole chunks
    for (i = 0; i < n - 4; i += 4) {
        packed[j++] = kjg_2bit_pack_unit(&unpacked[i]);
    }

    // pack the last chunk
    uint8_t remainder[4] = { 0, 0, 0, 0 };
    for (; i < n; i++) {
        remainder[i % 4] = unpacked[i];
    }
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

inline size_t kjg_2bit_packed_tda (const size_t n) {
    return ((n + 3) / 4);
}
