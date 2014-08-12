#include <stdint.h>
#include <string.h>

#define U1(n)     n,            n+1,              n+2,              n+3
#define U2(n) U1(n),      U1(n+256),      U1(n+256*2),      U1(n+256*3)
#define U3(n) U2(n),    U2(n+65536),    U2(n+65536*2),    U2(n+65536*3)
static const uint32_t UNPACK_LOOKUP[256] = {
        U3(0), U3(16777216), U3(16777216*2), U3(16777216*3)
};

#define P0(a, b, c, d) ((a*4 + b)*4 + c)*4 + d
#define P1(a, b, c)    { P0(a,b,c,0), P0(a,b,c,1), P0(a,b,c,2), P0(a,b,c,3) }
#define P2(a, b)       {   P1(a,b,0),   P1(a,b,1),   P1(a,b,2),   P1(a,b,3) }
#define P3(a)          {     P2(a,0),     P2(a,1),     P2(a,2),     P2(a,3) }
static const uint8_t PACK_LOOKUP[4][4][4][4] = { P3(0), P3(1), P3(2), P3(3) };

void kjg_2bit_pack(const size_t n, const uint8_t* unpacked, uint8_t* packed);

void kjg_2bit_unpack(const size_t n, const uint8_t* packed, uint8_t* unpacked);
