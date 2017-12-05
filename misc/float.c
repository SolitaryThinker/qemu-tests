#include <stdio.h>
#include <math.h>

/* little endian since we are emulating x86 */
#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x



static const double
ln2_hi  =  6.93147180369123816490e-01,/* 3fe62e42 fee00000 */
ln2_lo  =  1.90821492927058770002e-10,/* 3dea39ef 35793c76 */
two54   =  1.80143985094819840000e+16,  /* 43500000 00000000 */
Lp1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lp2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lp3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lp4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lp5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lp6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lp7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

static double zero = 0.0;

double log1p_l(double x)
{
  double hfsq,f,c,s,z,R,u;
  int k,hx,hu,ax;

  hx = __HI(x);/* high word of x */
    ax = hx&0x7fffffff;

  k = 1;
  if (hx < 0x3FDA827A) {/* x < 0.41422  */
    if(ax>=0x3ff00000) {/* x <= -1.0 */
      if(x==-1.0) return -two54/zero; /* log1p(-1)=+inf */
      else return (x-x)/(x-x);/* log1p(x<-1)=NaN */
    }
    if(ax<0x3e200000) {/* |x| < 2**-29 */
      if(two54+x>zero/* raise inexact */
          &&ax<0x3c900000)    /* |x| < 2**-54 */
        return x;
      else
        return x - x*x*0.5;
    }
    if(hx>0||hx<=((int)0xbfd2bec3)) {
      k=0;f=x;hu=1;}/* -0.2929<x<0.41422 */
  } 
  if (hx >= 0x7ff00000) return x+x;
  if(k!=0) {
    if(hx<0x43400000) {
      u  = 1.0+x; 
      hu = __HI(u);/* high word of u */
        k  = (hu>>20)-1023;
      c  = (k>0)? 1.0-(u-x):x-(u-1.0);/* correction term */
      c /= u;
    } else {
      u  = x;
      hu = __HI(u);/* high word of u */
        k  = (hu>>20)-1023;
      c  = 0;
    }
    hu &= 0x000fffff;
    if(hu<0x6a09e) {
      __HI(u) = hu|0x3ff00000;/* normalize u */
    } else {
      k += 1; 
      __HI(u) = hu|0x3fe00000;/* normalize u/2 */
        hu = (0x00100000-hu)>>2;
    }
    f = u-1.0;
  }
  hfsq=0.5*f*f;
  if(hu==0) {/* |f| < 2**-20 */
    if(f==zero) if(k==0) return zero;  
    else {c += k*ln2_lo; return k*ln2_hi+c;}
    R = hfsq*(1.0-0.66666666666666666*f);
    if(k==0) return f-R; else
      return k*ln2_hi-((R-(k*ln2_lo+c))-f);
  }
  s = f/(2.0+f); 
  z = s*s;
  R = z*(Lp1+z*(Lp2+z*(Lp3+z*(Lp4+z*(Lp5+z*(Lp6+z*Lp7))))));
  if(k==0) return f-(hfsq-s*(hfsq+R)); else
    return k*ln2_hi-((hfsq-(s*(hfsq+R)+(k*ln2_lo+c)))-f);
}





static const double
Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */


double __ieee754_log(double x)
{
  double hfsq,f,s,z,R,w,t1,t2,dk;
  int k,hx,i,j;
  unsigned lx;

  hx = __HI(x);/* high word of x */
    lx = __LO(x);/* low  word of x */

    k=0;
  if (hx < 0x00100000) {/* x < 2**-1022  */
    if (((hx&0x7fffffff)|lx)==0) 
      return -two54/zero;/* log(+-0)=-inf */
        if (hx<0) return (x-x)/zero;/* log(-#) = NaN */
          k -= 54; x *= two54; /* subnormal number, scale up x */
    hx = __HI(x);/* high word of x */
  } 
  if (hx >= 0x7ff00000) return x+x;
  k += (hx>>20)-1023;
  hx &= 0x000fffff;
  i = (hx+0x95f64)&0x100000;
  __HI(x) = hx|(i^0x3ff00000);/* normalize x or x/2 */
    k += (i>>20);
  f = x-1.0;
  if((0x000fffff&(2+hx))<3) {/* |f| < 2**-20 */
    if(f==zero) if(k==0) return zero;  else {dk=(double)k;
      return dk*ln2_hi+dk*ln2_lo;}
    R = f*f*(0.5-0.33333333333333333*f);
    if(k==0) return f-R; else {dk=(double)k;
      return dk*ln2_hi-((R-dk*ln2_lo)-f);}
  }
  s = f/(2.0+f); 
  dk = (double)k;
  z = s*s;
  i = hx-0x6147a;
  w = z*z;
  j = 0x6b851-hx;
  t1= w*(Lg2+w*(Lg4+w*Lg6)); 
  t2= z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7))); 
  i |= j;
  R = t2+t1;
  if(i>0) {
    hfsq=0.5*f*f;
    if(k==0) return f-(hfsq-s*(hfsq+R)); else
      return dk*ln2_hi-((hfsq-(s*(hfsq+R)+dk*ln2_lo))-f);
  } else {
    if(k==0) return f-s*(f-R); else
      return dk*ln2_hi-((s*(f-R)-dk*ln2_lo)-f);
  }
}





/* asinh(x)
 *  * Method :
 *  *MethodBased on 
 *  *MethodBasedasinh(x) = sign(x) * log [ |x| + sqrt(x*x+1) ]
 *  *sqrtwe have
 *  *haveasinh(x) := x  if  1+x*x=1,
 *  *haveasinh := sign(x)*(log(x)+ln2)) for large |x|, else
 *  *large := sign(x)*log(2|x|+1/(|x|+sqrt(x*x+1))) if|x|>2, else
 *         *if := sign(x)*log1p(|x| + x^2/(1 + sqrt(1+x^2)))  
 *          */


static const double 
one =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
ln2 =  6.93147180559945286227e-01, /* 0x3FE62E42, 0xFEFA39EF */
huge=  1.00000000000000000000e+300; 

double asinh_l(double x)
{
  double t,w;
  int hx,ix;
  hx = __HI(x);
  ix = hx&0x7fffffff;
  if(ix>=0x7ff00000) return x+x;/* x is inf or NaN */
    if(ix< 0x3e300000) {/* |x|<2**-28 */
      if(huge+x>one) return x;/* return x inexact except 0 */
    } 
  if(ix>0x41b00000) {/* |x| > 2**28 */
    w = __ieee754_log(fabs(x))+ln2;
  } else if (ix>0x40000000) {/* 2**28 > |x| > 2.0 */
    t = fabs(x);
    w = __ieee754_log(2.0*t+one/(sqrt(x*x+one)+t));
  } else {/* 2.0 > |x| > 2**-28 */
    t = x*x;
    w =log1p_l(fabs(x)+t/(one+sqrt(one+t)));
  }
  if(hx>0) return w; else return -w;
}

int main(){ /* compile with gcc -lm test.c */
  long double x, y;
  x = -2.1073424255447017e-08;
  // this fails
  //y = asinh(x);

  // this passes
  y = asinh_l(x);
  printf("y = %20.16Le\n",y);
  printf("should be -2.1073424255447017e-08\n");
}
