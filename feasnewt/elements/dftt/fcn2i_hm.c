#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.33333333333333333333333333333e-00        /* -4.0/3.0       */

#define d	1.0e-4					   /* delta          */
#define fd2     4.0e-8                                     /*  4.0*delta^2   */

int o_fcn_hm(double *obj, const double x[6])
{
  static double matr[4], f, t1, t2;
  static double fmat[3], g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];

  matr[2] = x[4] - x[3];
  matr[3] = x[5] - x[3];

  /* Calculate det(M). */
  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] - 1.0;
  fmat[1] = matr[0]*matr[2] + matr[1]*matr[3];
  fmat[2] = matr[2]*matr[2] + matr[3]*matr[3] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] +
                            fmat[2]*fmat[2];

  /* Calculate objective function. */
  (*obj) = 1.0 / (a * f * pow(g, b));
  return 0;
}

