#include <stdio.h>
#include <math.h>

int main(){ /* compile with gcc -lm test.c */
  long double x, y;
  x = -2.1073424255447017e-08;
  y = asinhl(x);
  printf("y = %20.16Le\n",y);
  printf("should be -2.1073424255447017e-08\n");
}
