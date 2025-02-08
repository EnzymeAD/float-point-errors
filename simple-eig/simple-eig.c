#include <fenv.h>
#include <math.h>
#include <stdint.h>
#define TRUE 1
#define FALSE 0

// ## PRE J2: 0, 0.33
// ## PRE J3: -0.002, 0.002
__attribute__((noinline))
double example(double J2, double J3) {
  // angle used to find eigenvalues
  double tmp = (0.5 * J3) * pow(3.0 / J2, 1.5);
  double alpha = acos(fmin(fmax(tmp, -1.0), 1.0)) / 3.0;

  // consider the most distinct eigenvalue first
  if (6.0 * alpha < M_PI) {
    return 2.0 * sqrt(J2 / 3.0) * cos(alpha);
  } else {
    return 2.0 * sqrt(J2 / 3.0) * cos(alpha + 2.0 * M_PI / 3.0);
  }
}