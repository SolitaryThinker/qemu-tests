// https://stackoverflow.com/questions/13725802/properties-of-80-bit-extended-precision-computations-starting-from-double-precis
//
#include <stdio.h>

double interpol_64(double u1, double u2, double u3)
{ 
  return u2 * (1.0 - u1) + u1 * u3;  
}

double interpol_80(double u1, double u2, double u3)
{ 
  return u2 * (1.0 - (long double)u1) + u1 * (long double)u3;  
}

int main()
{
  double y64,y80,u1,u2,u3;
  u1 = 0.025;
  u2 = 0.195;
  u3 = 0.195;
  y64 = interpol_64(u1, u2, u3);
  y80 = interpol_80(u1, u2, u3);
  printf("u2: %a\ny64:%a\ny80:%a\n", u2, y64, y80);
}
