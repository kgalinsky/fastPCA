#include "kjg_2bit.h"

#include <string.h>

// macros for unpacking array generation
#define U1(n)     n,            n+1,              n+2,              n+3
#define U2(n) U1(n),      U1(n+256),      U1(n+256*2),      U1(n+256*3)
#define U3(n) U2(n),    U2(n+65536),    U2(n+65536*2),    U2(n+65536*3)

// unpack lookup array
static const uint32_t UNPACK_LOOKUP[256] = { U3(0), U3(16777216), U3(
        16777216 * 2), U3(16777216 * 3) };

// packs 4 integers into a byte
#define PACK(a, b, c, d) ((a*4 + b)*4 + c)*4 + d

// macros for pack array generation
#define P1(a, b, c)    { PACK(a,b,c,0), PACK(a,b,c,1), PACK(a,b,c,2), PACK(a,b,c,3) }
#define P2(a, b)       {     P1(a,b,0),     P1(a,b,1),     P1(a,b,2),     P1(a,b,3) }
#define P3(a)          {       P2(a,0),       P2(a,1),       P2(a,2),       P2(a,3) }

// pack lookup array (faster than running PACK)
static const uint8_t PACK_LOOKUP[4][4][4][4] = { P3(0), P3(1), P3(2), P3(3) };

size_t kjg_2bit_unpack (
        const size_t n,
        const uint8_t* packed,
        uint8_t* unpacked) {
    size_t i, j = 0;
    for (i = 0; i < n - 4; i += 4) {
        memcpy(&unpacked[i], &UNPACK_LOOKUP[packed[j++]], 4);
    }
    memcpy(&unpacked[i], &UNPACK_LOOKUP[packed[j]], n - i);
    return (j);
}

size_t kjg_2bit_unpack_or (
        const size_t n,
        const uint8_t* packed,
        const uint8_t* mask,
        uint8_t* unpacked) {
    size_t i, j = 0;
    for (i = 0; i < n; i += 4) {
        memcpy(&unpacked[i], &UNPACK_LOOKUP[packed[j] | mask[j++]],
                sizeof(uint32_t));
    }
    return (j);
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

inline uint8_t kjg_2bit_pack_unit (const uint8_t* unpacked) {
    return (PACK_LOOKUP[unpacked[3]][unpacked[2]][unpacked[1]][unpacked[0]]);
}

inline size_t kjg_2bit_packed_tda (const size_t n) {
    return ((n + 3) / 4);
}
