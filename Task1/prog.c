#include <math.h>
#include <stdio.h>

#define C 4.29e-13

double phi (double x, int alpha) {
  return log( cosh(((1.) / alpha) + ((1.) / ((1 + cos (x) * cos(x))  ))  )) - log(alpha);
}

int main(void) {
  double x = 0;
  double y = 0;
  for (int i = 1; i <= 10; i++) {
    y = phi(x, i);
    while (fabs(x - y) > C) {
      x = phi(x, i);
      y = phi(y, i);
    }
    printf("alpha = %d, x =  %lf\n", i, y);
    
  }
  return 0;
}
