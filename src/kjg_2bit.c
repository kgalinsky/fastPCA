#include "kjg_2bit.h"

void kjg_2bit_pack(const size_t n, const uint8_t* unpacked, uint8_t* packed) {
    size_t i, j = 0;
    for (i = 3; i < n; i += 7) {
        packed[j++] =
                PACK_LOOKUP[unpacked[i--]][unpacked[i--]][unpacked[i--]][unpacked[i]];
    }
}

void kjg_2bit_unpack(const size_t n, const uint8_t* packed, uint8_t* unpacked) {
    size_t i;
    for (i = 0; i < n; i += 4) {
        memcpy(&unpacked[i], &UNPACK_LOOKUP[packed[i / 4]], sizeof(uint32_t));
    }
}
