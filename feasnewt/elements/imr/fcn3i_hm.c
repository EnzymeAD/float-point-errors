#include <math.h>
#include "fcn.h"

#define a       4.16666666666666666666666666666e-02        /*  1.0/24.0     */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/ 3.0      */

int o_fcn_hm(double *obj, const double x[12])
{
  static double matr[9], f;
  static double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

  /* Calculate det(M). */

  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  (*obj) = 1.0 / (a * f * pow(g, b));
  return 0;
}

