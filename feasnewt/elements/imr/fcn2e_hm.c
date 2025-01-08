#include <math.h>
#include "fcn.h"

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */

int o_fcn_hm(double *obj, const double x[6])
{
  static double matr[4], f;
  static double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[4] - x[3];
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  /* Calculate objective function. */
  (*obj) = 1.0 / (a * f * pow(g, b));
  return 0;
}

