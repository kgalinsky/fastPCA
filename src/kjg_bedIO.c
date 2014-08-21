#include "kjg_bedIO.h"

// macros for packing array generation
#define P1(a, b, c) { KJG_BEDIO_PACK(a,b,c,0), \
                      KJG_BEDIO_PACK(a,b,c,1), \
                      KJG_BEDIO_PACK(a,b,c,2), \
                      KJG_BEDIO_PACK(a,b,c,3) }

#define P2(a, b)    { P1(a,b,0), P1(a,b,1), P1(a,b,2), P1(a,b,3) }
#define P3(a)       { P2(a,0),   P2(a,1),   P2(a,2),   P2(a,3) }
#define P4

const uint8_t KJG_BEDIO_PACK_LOOKUP[4][4][4][4] = { P3(0), P3(1), P3(2), P3(3) };

// macros for unpacking array generation
#define U1(n) KJG_BEDIO_UNPACK(n), \
              KJG_BEDIO_UNPACK(n | 1), \
              KJG_BEDIO_UNPACK(n | 2), \
              KJG_BEDIO_UNPACK(n | 3)

#define U2(n) U1(n), U1(n | (1 << 2)), U1(n | (2 << 2)), U1(n | (3 << 2))
#define U3(n) U2(n), U2(n | (1 << 4)), U2(n | (2 << 4)), U2(n | (3 << 4))

const uint8_t KJG_BEDIO_UNPACK_LOOKUP[256][4] = { U3(0), U3((1<<6)), U3((2 << 6)), U3((3 << 6)) };

// macros for geno2bed array generation
#define G1(n) KJG_BEDIO_GENO2BED((n)), \
              KJG_BEDIO_GENO2BED((n | 1)), \
              KJG_BEDIO_GENO2BED((n | 2)), \
              KJG_BEDIO_GENO2BED((n | 3))

#define G2(n) G1(n), G1(n | (1 << 2)), G1(n | (2 << 2)), G1(n | (3 << 2))
#define G3(n) G2(n), G2(n | (1 << 4)), G2(n | (2 << 4)), G2(n | (3 << 4))

const uint8_t KJG_BEDIO_GENO2BED_LOOKUP[256] = { G3(0), G3((1 << 6)), G3((2 << 6)), G3((3 << 6)) };

// macros for geno2bed array generation
#define B1(n) KJG_BEDIO_BED2GENO(n), \
              KJG_BEDIO_BED2GENO(n | 1), \
              KJG_BEDIO_BED2GENO(n | 2), \
              KJG_BEDIO_BED2GENO(n | 3)

#define B2(n) B1(n), B1(n | (1 << 2)), B1(n | (2 << 2)), B1(n | (3 << 2))
#define B3(n) B2(n), B2(n | (1 << 4)), B2(n | (2 << 4)), B2(n | (3 << 4))

const uint8_t KJG_BEDIO_BED2GENO_LOOKUP[256] = { B3(0), B3((1 << 6)), B3((2 << 6)), B3((3 << 6)) };
