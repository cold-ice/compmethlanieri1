#include <stdint.h>

//Q15 multiply
// x and y are 16-bit integers (Q15)
//#define FIX_MPY(x,y)  ((int16_t)(((int32_t)x * (int32_t)y)>>15))
int16_t FIX_MPY(int16_t x, int16_t y) {

	return ((int16_t)(((int32_t)x * (int32_t)y)>>15));
}

//25-bit by 18-bit multiply, with output scaled to 25-bits, Q24xQ17
// x is 25-bits stored as 32, y is 18-bits stored as 32
int32_t FIX_MPY25by18(int32_t x,int32_t y) {

  return ((int32_t)(((long long)x * (long long)y)>>17));
}

//Saturated addition for Q15
int16_t SAT_ADD16(int16_t x,int16_t y) {

  if ((int32_t)x + (int32_t)y > 32767)
    return(32767);
  else if ((int32_t)x + (int32_t)y < -32767)
    return(-32768);
  else
    return(x+y);
}

//Saturated addition for Q24
int32_t SAT_ADD25(int32_t x,int32_t y) {

  if (x + y > (1<<24)-1)
    return((1<<24)-1);
  else if (x + y < -(1<<24))
    return(-(1<<24));
  else
    return(x+y);
}
