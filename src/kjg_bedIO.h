#ifndef KJG_BEDIO_H_
#define KJG_BEDIO_H_

#include <stdio.h>

#include "kjg_2bit.h"

typedef struct {
    const size_t m;   // number of SNPs
    const size_t n;   // number of samples
    FILE* stream;
} kjg_bedIO;

// unpacked integer conversion tables
#define KJG_BEDIO_123_231(n) ( (n) == 0 ? 0 : (n) == 3 ? 1 : ((n) + 1) )
#define KJG_BEDIO_231_123(n) ( (n) == 0 ? 0 : (n) == 1 ? 3 : ((n) - 1) )

// pack macros
#define KJG_BEDIO_PACK(a, b, c, d) \
    (  KJG_BEDIO_123_231(a) | \
      (KJG_BEDIO_123_231(b) << 2) | \
      (KJG_BEDIO_123_231(c) << 4) | \
      (KJG_BEDIO_123_231(d) << 6) )

// unpack macros
#define KJG_BEDIO_UA(n) KJG_BEDIO_231_123( (n)       & 3)
#define KJG_BEDIO_UB(n) KJG_BEDIO_231_123(((n) >> 2) & 3)
#define KJG_BEDIO_UC(n) KJG_BEDIO_231_123(((n) >> 4) & 3)
#define KJG_BEDIO_UD(n) KJG_BEDIO_231_123(((n) >> 6) & 3)

#define KJG_BEDIO_UNPACK(n) { KJG_BEDIO_UA(n), \
                              KJG_BEDIO_UB(n), \
                              KJG_BEDIO_UC(n), \
                              KJG_BEDIO_UD(n) }

// conversion macros
#define KJG_BEDIO_GENO2BED(n) KJG_BEDIO_PACK( KJG_2BIT_UA(n), \
                                              KJG_2BIT_UB(n), \
                                              KJG_2BIT_UC(n), \
                                              KJG_2BIT_UD(n) )

#define KJG_BEDIO_BED2GENO(n) KJG_2BIT_PACK( KJG_BEDIO_UA(n), \
                                             KJG_BEDIO_UB(n), \
                                             KJG_BEDIO_UC(n), \
                                             KJG_BEDIO_UD(n) )

// lookup tables
extern const uint8_t KJG_BEDIO_PACK_LOOKUP[4][4][4][4];
extern const uint8_t KJG_BEDIO_UNPACK_LOOKUP[256][4];
extern const uint8_t KJG_BEDIO_GENO2BED_LOOKUP[256];
extern const uint8_t KJG_BEDIO_BED2GENO_LOOKUP[256];

#endif /* KJG_BEDIO_H_ */
