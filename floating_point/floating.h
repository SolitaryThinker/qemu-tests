#include <stdint.h>

#ifndef __FLOATING_H
#define __FLOATING_H

typedef uint64_t flag;
typedef uint64_t bits64;

typedef struct {
    uint64_t low;
    uint16_t high;
} floatx80;

/* extended precision helpers */

flag get_floatx80_sign(floatx80 a);
int32_t get_floatx80_exp(floatx80 a);
bits64 get_floatx80_sig(floatx80 a);

int32_t compact_floatx80(int32_t exp, bits64 sig);

#endif
