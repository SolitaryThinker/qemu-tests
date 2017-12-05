#include "floating.h"


flag get_floatx80_sign(floatx80 a) {
    return  a.high & 0x8000;
}

int32_t get_floatx80_exp(floatx80 a) {
    return a.high & 0x7FFF;
}

bits64 get_floatx80_sig(floatx80 a) {
    return a.low;
}

int32_t compact_floatx80(int32_t exp, bits64 sig) {
    return exp << 16 | sig >> 48;
}

//int main() {
    //int i = 1;
    //long double d = 1.0;
    //double d2 = 1.0;
    //float f = 1.0;

    //return 0;
//}
